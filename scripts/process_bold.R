#!/usr/bin/env Rscript
smth<-1
nrandrun<-1
ndvis<-8
zval<-( -1 )
clustval<-20
timefilter<-FALSE
residdesign<-FALSE
CORAC<-FALSE
dounivar<-FALSE 
sparval <- 0.1
sparval2<- -0.3
domotor <- FALSE
library(getopt)
initlist<-list()
options(digits=3)
Args <- commandArgs()
self<-Args[4]
self<-substring(self,8,nchar(as.character(self)))
spec = c(
'subjectid', 's', 1, "character" ," unique ID for subject ",
'design', 'd', 1, "character" ," design matrix ",
'bold', 'b', 1, "character" ," BOLD image ",
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
if ( is.null( opt$bold ) ) {
  fmri<-antsImageRead('bold.nii.gz',4)
} else {
  fmri<-antsImageRead( opt$bold ,4)
}
# fmrit<-antsImageClone(fmri)
ImageMath(4,fmri,'SliceTimingCorrection',fmri,0)
if ( opt$design == "task003" ) domotor<-TRUE
print(paste("Do Motor",domotor,opt$design))
tr<-as.numeric( opt$tr )
nbadframes<-5
blockfing = c(0, 36, 72, 108,144)+nbadframes
blockfoot <- blockfing + 12
blockmout <- blockfoot + 12 
blocko = c(1, 24, 48, 72, 96, 120, 144 )+nbadframes # covert 
if ( domotor ) ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfing , durations=rep(  6,  length( blockfing ) ) ,  rt=tr ) else ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blocko , durations=rep(  12,  length( blocko ) ) ,  rt=tr ) 
ohrf2 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfoot , durations=rep(  6,  length( blockfoot ) ) ,  rt=tr )
ohrf3 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockmout , durations=rep(  6,  length( blockmout ) ) ,  rt=tr )
if ( domotor ) ohrf<-cbind( ohrf , ohrf2, ohrf3 ) else ohrf<-cbind( ohrf )
ohrf[1:4,]<-0 # first few frames are junk
write.csv(ohrf,'antsr_hrf.csv',quote=F,row.names=F)
if ( ! file.exists("mat.mha") |  ! file.exists("mask.nii.gz") | ! file.exists("nuis.csv")  ) {
  myvarsin<-getfMRInuisanceVariables( fmri, moreaccurate = TRUE ,  maskThresh=100 )
  mat <- myvarsin$matrixTimeSeries
  antsImageWrite( as.antsImage(mat) , "motion_corrected.nii.gz" )
  antsImageWrite( as.antsImage(mat) , "mat.mha" )
  antsImageWrite( myvarsin$mask , "mask.nii.gz" )
  antsImageWrite( myvarsin$avg ,  paste("avg.nii.gz",sep="" )  )
  write.csv( data.frame( myvarsin$globalsignal, myvarsin$nuisancevariables ),  'nuis.csv' )
} else {
  avg <- new("antsImage", dimension = 3, pixeltype = fmri@pixeltype)
  antsMotionCorr(list(d = 3, a = fmri, o = avg))
  avg<-antsImageWrite( avg, 'avg.nii.gz')
  nuis<-read.csv('nuis.csv')
  mat<-as.matrix( antsImageRead( 'mat.mha', 2 ) )
  mask<-antsImageRead( 'mask.nii.gz', 3 )
  avg<-antsImageRead( 'avg.nii.gz', 3 )
  myvarsin<-list(  matrixTimeSeries = mat , globalsignal = nuis$myvarsin.globalsignal, nuisancevariables = nuis, mask=mask )
}
if ( domotor ) myvarsin$nuisancevariables<-cbind( myvarsin$nuisancevariable, hrf2=ohrf[,2], hrf3=ohrf[,3] )
mask<-myvarsin$mask
randuvar<-antsImageClone( mask )
randuvar[ mask > 0 ]<-0
randsccan<-antsImageClone( mask )
randsccan[ mask > 0 ]<-0
for ( randrun in 1:nrandrun ) {
  myvars<-myvarsin
  hrf<-ohrf
  mat<-myvars$matrixTimeSeries
  mat<-mat[5:nrow(hrf),]
  myvars$globalsignal<-myvars$globalsignal[ 5:nrow(hrf) ]
  myvars$nuisancevariables<-myvars$nuisancevariables[5:nrow(hrf),]
  hrf<-as.matrix( hrf[ 5:nrow(hrf), ] )
  sfmri<-antsImageClone( fmri )
  SmoothImage(4,fmri,3,sfmri)
  mat<-timeseries2matrix( sfmri, myvars$mask )
  mat<-whiten( mat[5:nrow(mat),] )
  mysubset<-c(1:(nrow(mat)*1))
  myruns<-c(rep(1,ndvis-1),2)
  myruns<-sample(myruns)
  myruns2<-c()
  nreps<-round( nrow(mat)/ndvis )-1
  for ( i in myruns ) myruns2<-c( myruns2, rep( i , nreps ) ) 
  if (  nrandrun > 1 ) mysubset<-which( myruns2 == 1 )
#  mysubset<-c(1:(nrow(mat)*1))
#  if ( nrandrun > 1 ) mysubset<-sample(mysubset)[1:round(length(mysubset)*0.8)]
  mat<-mat[mysubset,]
  myvars$nuisancevariables<-myvars$nuisancevariables[mysubset,]
  myvars$globalsignal<-myvars$globalsignal[mysubset]
  if ( domotor ) hrf<-hrf[mysubset,] else  hrf<-as.matrix(hrf[mysubset])
  if ( timefilter ) mat  <- filterfMRIforNetworkAnalysis( cbind(mat) , tr, cbfnetwork = "BOLD" , freqLo=0.01 , freqHi = 0.2  )$filteredTimeSeries
  if ( domotor ) colnames(hrf)<-c("hrf","hrf2","hrf3") else colnames(hrf)<-c("hrf")
  gsig<-myvars$globalsignal
  arval<-ar(gsig,FALSE,2)$ar
  arval<-shift(gsig,1)*arval[1]+shift(gsig,2)*arval[2]
  baseform<-" motion1 + motion2 + motion3 + compcorr1 + compcorr2 + compcorr3 + metricnuis"
  myform<-paste(baseform," +  globalsignal ")
  if ( CORAC ) myform<-paste(myform," + arval ")
  if ( domotor ) myform<-paste(myform," + hrf2 + hrf3 ")
  print( myform )
  if ( dounivar ) {
    fmrimodel <- taskFMRI( mat , hrf, myvars  , correctautocorr = FALSE,
                          residualizedesignmatrix  = residdesign, myformula=myform ) # 
    betas<-fmrimodel$beta
    maxbeta<-max( betas )
    betaimg<-antsImageClone( myvars$mask ) # put beta vals in image space
    betaimg[ myvars$mask > 0.5 ] <- betas
    betathresh<-min(betas[rev(order(betas))][1:( round( length(betas) * sparval ))])
                                        # ( mean( betas[betas>0.1] )+ 2.0 * sd(  betas[betas>0.1] ) )
    antsImageWrite(betaimg, paste("betas.nii.gz",sep="") )
    betasthresh<-betas
    betasthresh[ betas < betathresh ]<-0
    randuvar[ myvars$mask > 0.5 ] <-randuvar[ myvars$mask > 0.5 ] + betasthresh/nrandrun
    antsImageWrite(randuvar, paste("betast.nii.gz",sep="") )
    amat<-fmrimodel$fmrimat
############################################
    mydf<-data.frame(TrainMaxBeta=maxbeta)
    write.csv( mydf , file = "maxbeta.csv" )
############################################
  }# now some multivariate stuff
  mat<-( as.matrix( antsImageRead( 'mat.mha', 2 ) ) )
  hrf<-ohrf
  hrf<-hrf[5:nrow(mat),]
  mat<-mat[5:nrow(mat),]
  if ( timefilter ) mat  <- filterfMRIforNetworkAnalysis( cbind(mat) , tr, cbfnetwork = "BOLD" , freqLo=0.01 , freqHi = 0.2  )$filteredTimeSeries
  ndf<-data.frame( globalsignal= myvarsin$globalsignal,  myvarsin$nuisancevariables )
  mat<-mat[mysubset,]
  if ( domotor ) hrf<-hrf[mysubset,] else hrf<-as.matrix(hrf[mysubset])
  ndf<-ndf[mysubset,]
  if ( residdesign ) ndf<-data.frame( residuals( lm( as.matrix( ndf ) ~ hrf[,1] ) ) )
  mat<-residuals( lm( as.formula( paste(" mat ~ ",myform) ) , data = ndf ) )
  print(paste("done residualizing for sccan, randrun:",randrun))
  mypreds<-as.matrix( cbind( hrf, as.numeric(  hrf[,1] > 0 )  ) )
  mypreds<-as.matrix( cbind( hrf  ) )
#  if ( CORAC ) mat<-amat
  docca<-TRUE
  if ( docca ) {
  bestp<-1
  nv<-ncol(mypreds)+1
  sccan<-sparseDecom2( inmatrix=list( mat , mypreds ), inmask = c( myvars$mask , NA ) ,
                    sparseness=c( sparval , sparval2 ), nvecs=nv, its=4, smooth=smth,
                      perms=0, cthresh = c(clustval, 0) , robust=0, mycoption=0,
                      z=zval, initializationList = initlist )
  if ( length( initlist ) < nv ) { initlist<-sccan$eig1 }
  for (  kk in 1:nv ) {
    pv <-  cor.test( hrf[,1], sccan$projections[,kk] )$p.value
    if ( pv < bestp )
      {
      print(paste("choose",kk,pv))
      bestp<-pv
      myeig<-antsImageClone( sccan$eig1[[kk]] )
      }
    }
  print( sccan$eig2 )
  }
  if ( ! docca ) {
    bestp<-1
    nv<-20
#    if ( length( initlist ) > 0 ) sparval <- 0 
    sccan<-sparseDecom( inmatrix= mat , inmask = mask  , z=zval,
                    sparseness=c( sparval  ), nvecs=nv, its=3, smooth=1,
                       cthresh = c(10)  , initializationList = initlist )
    if ( length( initlist ) < nv ) initlist<-sccan$eigenanatomyimages
    for ( jj in 1:ncol(sccan$projections) )
      {
        vox<-sccan$projections[,jj]
        locpval<-cor.test( vox, hrf)$p.value
        if ( locpval <  bestp )
          {
            bestp<-locpval
            myeig<-antsImageClone( sccan$eigenanatomyimages[[jj]] )
            print( paste( jj , cor.test( vox, hrf)$p.value ) )
          }
#    myuform<-as.formula( paste( " vox ~ hrf + " , myform ) )
#    print( summary( lm( vox ~ hrf , data = ndf ) ) )
      }
  }
  ImageMath(3,myeig,'abs',myeig)
  ImageMath(3,randsccan,'+',randsccan,myeig)
  antsImageWrite( randsccan ,  paste(opt$output,"sccan.nii.gz",sep="")  )
}
randmean<-mean( randsccan[ randsccan > 1.e-8 ] )
randsd<-sd( randsccan[ randsccan > 1.e-8 ] )
esparval<-sum( randsccan > (randmean+randsd) ) / ncol(mat) 
print( paste( "estimated sparseness:", esparval , " with thresh " , (randmean+randsd)) )
#############################################################################
q()
myvars<-myvarsin
hrf<-ohrf
mat<-myvars$matrixTimeSeries
mat<-mat[5:nrow(hrf),]
myvars$globalsignal<-myvars$globalsignal[ 5:nrow(hrf) ]
myvars$nuisancevariables<-myvars$nuisancevariables[5:nrow(hrf),]
if ( timefilter ) mat  <- filterfMRIforNetworkAnalysis( cbind(mat) , tr, cbfnetwork = "BOLD" , freqLo=0.01 , freqHi = 0.2  )$filteredTimeSeries
ndf<-data.frame( globalsignal = myvars$globalsignal,  myvars$nuisancevariables )
if ( residdesign ) ndf<-data.frame( residuals( lm( as.matrix( ndf ) ~ hrf[,1] ) ) )
mat<-residuals( lm( as.formula( paste(" mat ~ ",myform) ) , data = ndf ) )
cblock <- as.numeric( hrf ) 
mypreds<-as.matrix( cbind( cblock, as.numeric( cblock > 0 )  ) )
sccan<-sparseDecom2( inmatrix=list( mat , mypreds ), inmask = c( myvars$mask , NA ) ,
                    sparseness=c( esparval*(1) , -1 ), nvecs=2, its=5, smooth=1,
                      perms=0, cthresh = c(0, 0) , robust=0 ) 
ImageMath(3,sccan$eig1[[1]],'abs',sccan$eig1[[1]])
antsImageWrite( sccan$eig1[[1]] ,  paste("sccan.nii.gz",sep="")  )
