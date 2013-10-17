#!/usr/bin/env Rscript
library(getopt)
options(digits=3)
Args <- commandArgs()
self<-Args[4]
self<-substring(self,8,nchar(as.character(self)))
spec = c(
'boldimg', 'b', 1, "character" ," testing bold matrix ",
'activation', 'a', 1, "character" ," training 3D activation image ",
'hrf', 'h', 1, "character" ," hrf csv file ",
'mask', 'm', 1, "character" ," 3D mask antsImage ",
'output', 'o', 1, "character" ," output ")
spec=matrix(spec,ncol=5,byrow=TRUE)
opt = getopt(spec)
if ( is.null( opt$boldimg ) ) {
  print( 'no 2D bold matrix , quitting , use --boldimg option')
q()
}
if ( is.null( opt$hrf ) ) {
  print( 'no hrf , quitting , use --hrf option')
q()
}
if ( is.null( opt$mask ) ) {
  print( 'no 3D bold mask , quitting , use --mask option')
q()
}
if ( is.null( opt$activation ) ) {
  print( 'no test activation (3D) image , quitting , use --activation option')
q()
}
if ( is.null( opt$output ) ) {
  print( 'no output file name , quitting , use --output option')
q()
}
######
library(ANTsR)
acti<-antsImageRead(opt$activation,3)
mat<-as.matrix( antsImageRead(opt$boldimg,2) )
if ( ! file.exists( opt$hrf ) )
  {
  print( paste(opt$hrf," file does not exist " ) )
  q()
  }
hrf <- as.numeric( read.csv( opt$hrf )$V1 )
mat<-mat[5:length(hrf),]
hrf<-hrf[ 5:length(hrf) ]

if ( FALSE ) {
  fmri<-antsImageRead(opt$boldimg,4)
  ImageMath(4,fmri,'SliceTimingCorrection',fmri,'bspline')
  myvars<-getfMRInuisanceVariables( fmri, moreaccurate = FALSE ,  maskThresh=100 )
  mat <- myvars$matrixTimeSeries
  mat<-mat[5:length(hrf),]
  myvars$globalsignal<-myvars$globalsignal[ 5:length(hrf) ]
  myvars$nuisancevariables<-myvars$nuisancevariables[5:length(hrf),]
  hrf<-hrf[ 5:length(hrf) ]
# mat  <- filterfMRIforNetworkAnalysis( cbind(mat) , tr, cbfnetwork = "BOLD" , freqLo=0.001 , freqHi = 0.15  )$filteredTimeSeries
  myform<-"motion1 + motion2 + motion3 + compcorr1 + compcorr2 + compcorr3 + globalsignal "
# now some multivariate stuff 
  ndf<-data.frame( globalsignal = myvars$globalsignal,  myvars$nuisancevariables )
  ndf<-data.frame( residuals( lm( as.matrix( ndf ) ~ hrf ) ) )
  mat<-residuals( lm( as.formula( paste(" mat ~ ",myform) ) , data = ndf ) )
  print("done residualizing, now testing")
}
###########################################################################
mask<-antsImageRead(opt$mask,3)
acti<-abs( acti[ mask > 0.5 ] )
mymean<-mean( acti[ acti > 1.e-10 ] )
mysd<-sd( acti[ acti > 1.e-10 ]  )
thresh<-(mymean+mysd*1)
ww<-which( acti >  thresh )
proj<-mat[,ww] %*% acti[ww]
print(paste(" Correlation of projection ",length(ww)," at thresh ", thresh ))
mycorr<-cor.test( hrf, proj )
print( mycorr )
mycorrs <- rep( 1, sum( mask > 0.5) )
myqvals <- rep( 1, sum( mask > 0.5) )
for ( nc in ww ) {
  mycorrs[nc]<-cor.test( hrf, mat[,nc] )$p.value
}
myqvals[ ww ]<-p.adjust( mycorrs[ ww ] , method="BH" )
print(paste( " Over all ", length(ww) , " voxels ") )
print( paste( mean(myqvals[ ww ]), min( myqvals[ ww ] ), max( myqvals[ ww ] ) ) )
outputdf<-data.frame(ProjectionPVal=mycorr$p.value,VoxMeanPval=mean(myqvals[ ww ]),VoxMinPval=min( myqvals[ ww ] ),VoxMaxPval=max( myqvals[ ww ] ), nvox=length(ww) )
write.csv(outputdf, paste( opt$output,'.csv' ,sep='' ) ,row.names=F,quote=F)
acti<-antsImageClone( mask )
acti[ mask > 0.5 ]<-0
acti[ mask > 0.5 ]<-1-myqvals
antsImageWrite( acti, paste( opt$output,'.nii.gz' ,sep='' ) )

