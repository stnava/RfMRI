#!/usr/bin/env Rscript
library(getopt)
options(digits=3)
Args <- commandArgs()
self<-Args[4]
self<-substring(self,8,nchar(as.character(self)))
spec = c(
'subjectid', 's', 1, "character" ," unique ID for subject ",
'design', 'd', 1, "character" ," design matrix ",
'tr', 't', 1, "character" ," BOLD tr ",
'run', 'n', 1, "character" ," which run ",
'output', 'o', 1, "character" ," output ")
spec=matrix(spec,ncol=5,byrow=TRUE)
opt = getopt(spec)
if ( is.null( opt$design ) ) {
  print( 'no design , quitting , use --design option')
q()
}
if ( is.null( opt$tr ) ) {
  print( 'no tr , quitting , use --tr option')
q()
}
if ( is.null( opt$run ) ) {
  print( 'no run , quitting , use --run option')
q()
}
if ( is.null( opt$subjectid ) ) {
  print( 'no subjectid , quitting , use --subjectid option')
q()
}

library(ANTsR)
sparval <- 0.05
fmri<-antsImageRead('bold.nii.gz',4)
ImageMath(4,fmri,'SliceTimingCorrection',fmri,'bspline')
stim<-as.numeric(read.table(opt$design)$V1)
tr<-as.numeric( opt$tr )
onsets<-round(stim/tr)
print( onsets )
blockfoot = c(12, 48, 84, 120, 156)
blockfing = c(0, 36, 72, 108,144)
hrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfing , durations=rep(  12,  length( blockfing ) ) ,  rt=2.5 )
hrf[1:4]<-0 # first few frames are junk 
myvars<-getfMRInuisanceVariables( fmri, moreaccurate = TRUE ,  maskThresh=100 )
mat <- myvars$matrixTimeSeries
mat<-mat[5:length(hrf),]
myvars$globalsignal<-myvars$globalsignal[ 5:length(hrf) ]
myvars$nuisancevariables<-myvars$nuisancevariables[5:length(hrf),]
hrf<-hrf[ 5:length(hrf) ]
mat  <- filterfMRIforNetworkAnalysis( cbind(mat) , tr, cbfnetwork = "BOLD" , freqLo=0.01 , freqHi = 0.2  )$filteredTimeSeries
myform<-"motion1 + motion2 + motion3 + compcorr1 + compcorr2 + compcorr3 + globalsignal "
if ( TRUE ) {
  sfmri<-antsImageClone( fmri )
  SmoothImage(4,fmri,4,sfmri)
  smat<-timeseries2matrix( sfmri, myvars$mask )
  smat<-smat[5:nrow(smat),]
  fmrimodel <- taskFMRI( smat , hrf, myvars  , correctautocorr = T,
   residualizedesignmatrix  = F, myformula=myform ) # 
  betas<-fmrimodel$beta
  maxbeta<-max( betas )
  betaimg<-antsImageClone( myvars$mask ) # put beta vals in image space
  betaimg[ myvars$mask > 0.5 ] <- betas
  betathresh<-min(betas[rev(order(betas))][1:( round( length(betas) * sparval ))])
    # ( mean( betas[betas>0.1] )+ 2.0 * sd(  betas[betas>0.1] ) )
  antsImageWrite(betaimg, paste("betas.nii.gz",sep="") )
  betasthresh<-betas
  betasthresh[ betas < betathresh ]<-0
  betaimg[ myvars$mask > 0.5 ] <-betasthresh
  antsImageWrite(betaimg, paste("betast.nii.gz",sep="") )
  amat<-fmrimodel$fmrimat
}
# now some multivariate stuff 
ndf<-data.frame( globalsignal = myvars$globalsignal,  myvars$nuisancevariables )
ndf<-data.frame( residuals( lm( as.matrix( ndf ) ~ hrf ) ) )
mat<-residuals( lm( as.formula( paste(" mat ~ ",myform) ) , data = ndf ) )
antsImageWrite( as.antsImage(mat) , "mat.mha" )
antsImageWrite( myvars$mask , "mask.nii.gz" )
write.csv(hrf,'antsr_hrf.csv',quote=F,row.names=F)
print("done residualizing for sccan")
cblock <- as.numeric(hrf) 
mypreds<-as.matrix( cbind( cblock, as.numeric( cblock > 0 )  ) )
sccan<-sparseDecom2( inmatrix=list( mat , mypreds ), inmask = c( myvars$mask , NA ) ,
  sparseness=c( sparval , 1 ), nvecs=ncol(mypreds), its=5, smooth=1,
  perms=2, cthresh = c(5, 0) , robust=0 ) 
antsImageWrite( myvars$avg ,  paste("avg.nii.gz",sep="" )  )
ImageMath(3,sccan$eig1[[1]],'abs',sccan$eig1[[1]])
ImageMath(3,sccan$eig1[[2]],'abs',sccan$eig1[[2]])
antsImageWrite( sccan$eig1[[1]] ,  paste("sccan.nii.gz",sep="" )  )
antsImageWrite( sccan$eig1[[2]] ,  paste("sccan2.nii.gz",sep="" )  )
############################################
outputstruct<-paste("TrainMaxBeta:",maxbeta)
cat( outputstruct , file = "maxbeta.txt" )
############################################
