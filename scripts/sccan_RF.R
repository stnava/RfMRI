library( ANTsR )
library( randomForest)
library( vegan )
library( mboost )
dowhiten <- TRUE
useglmb <- FALSE
takeoutresid <- TRUE
subnum<-"010"
print(paste("a random forest test on fmri ... takeoutresid? ",takeoutresid,'subject',subnum))
rootdir<-paste("/Users/stnava/data/data_gorgolewski/RfMRI/")
gimgfn<-paste(rootdir,"OptGroupTask003Run001sccan.nii.gz",sep='')
groupimg<-antsImageRead( gimgfn , 3 )
setwd(paste(rootdir,"/group_analysis/sub",subnum,"/task003/run001",sep=''))
fmri<-antsImageRead(paste('sub',subnum,'_group.nii.gz',sep=''),4)
tr<-2.5
nbadframes<-5
domotor<-TRUE
blockfing = c(0, 36, 72, 108,144)+nbadframes
blockfoot <- blockfing + 12
blockmout <- blockfoot + 12 
blocko = c(1, 24, 48, 72, 96, 120, 144 )+nbadframes # covert 
if ( domotor ) ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfing , durations=rep(  6,  length( blockfing ) ) ,  rt=tr ) else ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blocko , durations=rep(  12,  length( blocko ) ) ,  rt=tr ) 
ohrf2 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfoot , durations=rep(  6,  length( blockfoot ) ) ,  rt=tr )
ohrf3 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockmout , durations=rep(  6,  length( blockmout ) ) ,  rt=tr )
if ( domotor ) ohrf<-cbind( ohrf , ohrf2, ohrf3 ) else ohrf<-cbind( ohrf )
ohrf[1:4,]<-0 # first few frames are junk
temimg<-antsImageRead('../../../../template/template.nii.gz',3)
tem<-antsImageRead('../../../../template/aal.nii.gz',3)
mask<-antsImageClone( tem ) 
nuis<-read.csv('nuis.csv')
th<-c(1,90)
thresholdBool <- mask > th[2] | mask < th[1]
mask[   thresholdBool ]<-0
mask[ ! thresholdBool ]<-1
groupvec<-groupimg[ mask > 0 ]
mat<-timeseries2matrix( fmri , mask )
myglobsig<-apply( mat, FUN=mean, MARGIN=1 )
boldthresh<-1200 # mean(myglobsig)-6*sd(myglobsig)
selector<-myglobsig > boldthresh
selector[1:(nbadframes-1)] <- FALSE
mat<-subset( mat , selector )
nuis<-subset( nuis , selector )
myhrf<-subset( ohrf , selector )
nuis<-as.matrix( residuals( lm( as.matrix(nuis) ~ myhrf ) ) )
if ( takeoutresid )  mat<-( as.matrix( residuals( lm( mat  ~ as.matrix( nuis ) ) ) ) )
myglobsig<-apply( mat, FUN=mean, MARGIN=1 )
mat[,1]<-myglobsig
if ( dowhiten ) mat <- temporalwhiten( mat ) 
plot(   myhrf[,1], type='l' )
points( myhrf[,2], type='l', col='blue' )
points( myhrf[,3], type='l', col='red' )
matpreds<-matrix( rep(0,nrow(myhrf)*4 ) , nrow=nrow(myhrf) )
binpreds<-rep(0,nrow(myhrf))
mysums<-rep(0,nrow(myhrf))
for ( i in 5:nrow(myhrf) )
  {
  binpreds[i]<-which.max( myhrf[i,] )
  mysums[i]<-sum( abs( myhrf[i,] ) )
  if ( mysums[i] < 1.2   ) binpreds[i]<-0
  matpreds[i,binpreds[i]+1]<-1
  }
points( mysums, type='l', col='green' )
binpreds <- as.factor( binpreds )
myglobsigresid<-residuals( lm( myglobsig ~ myhrf ) )
threshinit<-0.3
img1<-antsImageClone( mask )
img1[ mask > 0 ]<-sparsify( myhrf[,1] %*% mat , threshinit, mask )
img2<-antsImageClone( mask )
img2[ mask > 0 ]<-sparsify( myhrf[,2] %*% mat , threshinit, mask )
img3<-antsImageClone( mask )
img3[ mask > 0 ]<-sparsify( myhrf[,3] %*% mat , threshinit, mask )
img4<-antsImageClone( mask )
img4[ mask > 0 ]<- sparsify( myglobsigresid %*% mat, threshinit, mask )
img5<-antsImageClone( mask )
img5[ mask > 0 ]<- sparsifyv( svd( mat )$u[2,] %*% mat, threshinit, mask )
initlist<-list( img1,img2,img3,img4 )
myhrfb<-matrix( as.numeric( myhrf > 0.5 ) , nrow=nrow(myhrf) )
mypreds <- cbind(myhrf,myhrfb)
ff<-sparseDecom2( inmatrix=list(mat,mypreds), inmask=list(mask,NA), perms=0, its=45, mycoption=1, sparseness=c( -0.02, 0.1 ) , nvecs= ncol( mypreds ), smooth=1, cthresh=c(10,0), ell1 = 11 , z=-1 ) #  , initializationList=initlist )
print( ff$eig2 )
mysccanimages<-imageListToMatrix( imageList=ff$eig1, mask=mask)
# ff<-sparseDecom( inmatrix=mat, inmask=mask,its=5,
#                  mycoption=1, sparseness=c(-0.1) ,nvecs=11, 
#                  smooth=1, cthresh=10, z=-1 )
# names(ff)[2]<-"eig1"
# names(ff)[3]<-"projections"
for ( i in 1:length(ff$eig1) ) {
  visimg<-antsImageClone( ff$eig1[[i]] )
  ImageMath(3,visimg,'abs',visimg)
  ImageMath(3,visimg,'Normalize',visimg)
  plotANTsImage( myantsimage=temimg, functional=list(visimg) , slices="2x25x1", ,axis=0, color=c("red") , threshold="0.01x1", outname=paste("temp",i,".jpg",sep='') ) 
}
gimgproj <- ( mat ) %*% groupvec
mysccanpreds <- ( mat  ) %*% t( mysccanimages )
colnames( mysccanpreds )<-colnames( ff$projections )
mydf<-data.frame( binpreds1 = myhrf[,1] , imgs = mysccanpreds  , gimgproj = gimgproj )
my.rf1<-randomForest( binpreds1 ~ . , data=mydf )
if ( useglmb ) my.rf1<-          glmboost( binpreds1 ~ . , data=mydf )
mydf<-data.frame( binpreds2 = myhrf[,2] , imgs = mysccanpreds  , gimgproj = gimgproj )
my.rf2<-randomForest( binpreds2 ~ . , data=mydf )
if ( useglmb ) my.rf2<-          glmboost( binpreds2 ~ . , data=mydf )
mydf<-data.frame( binpreds3 = myhrf[,3] , imgs = mysccanpreds  , gimgproj = gimgproj )
my.rf3<-randomForest( binpreds3 ~ . , data=mydf )
if ( useglmb ) my.rf3<-          glmboost( binpreds3 ~ . , data=mydf )
######################################################################
fmri2<-antsImageRead( paste('../run002/sub',subnum,'_group.nii.gz',sep='') , 4 )
nuis2<-read.csv('nuis.csv')
mat2<-timeseries2matrix( fmri2 , mask )
myglobsig2<-apply( mat2, FUN=mean, MARGIN=1 )
selector2<-myglobsig2 > boldthresh
selector2[1:(nbadframes-1)] <- FALSE
# selector2[50:length(selector2)] <- FALSE # for testing sub-sequence prediction
mat2<-subset( mat2 , selector2 )
nuis2<-subset( nuis2 , selector2 )
myhrf2<-subset( ohrf , selector2 )
nuis2<-as.matrix( residuals( lm( as.matrix(nuis2) ~ myhrf2 ) ) )
if ( takeoutresid ) mat2<-( as.matrix( residuals( lm( mat2  ~ as.matrix( nuis2 ) ) ) ) )
myglobsig2<-apply( mat2, FUN=mean, MARGIN=1 )
mat2[,1]<-myglobsig2
if ( dowhiten ) mat2 <- temporalwhiten( mat2 ) 
gimgproj <- ( mat2 )  %*% groupvec
mysccanpreds <- ( mat2  ) %*% t( mysccanimages )
colnames( mysccanpreds )<-colnames( ff$projections )
mydf2<-data.frame(  imgs = mysccanpreds , gimgproj = gimgproj )
########################################################################################################
mypred1<-predict( my.rf1 , newdata = mydf2 )
mypred2<-predict( my.rf2 , newdata = mydf2 )
mypred3<-predict( my.rf3 , newdata = mydf2 )
print( cor.test(  mypred1 , myhrf2[,1] ) )
plot(  myhrf2[,1] , type='l' )
points(   mypred1 , type='l' , col='green')
Sys.sleep(3)
print( cor.test(  mypred2 , myhrf2[,2] ) )
plot(  myhrf2[,2] , type='l' )
points(   mypred2 , type='l' , col='green')
Sys.sleep(3)
print( cor.test(  mypred3 , myhrf2[,3] ) )
plot(  myhrf2[,3] , type='l' )
points(   mypred3 , type='l' , col='green')
Sys.sleep(3)
pdf(paste(rootdir,'/figure/decoding',subnum,".pdf",sep=''),width=7,height=3)
plot(   mypred1 , type='l' , col='black',lty=2,ylim=c(-1,2),main="BOLD (dashed) decodes stimulus HRF (solid)")
points(   mypred2 , type='l' , col='red',lty=2)
points(   mypred3 , type='l' , col='green',lty=2)
points(   myhrf[,1] , type='l' , col='black',pch=2)
points(   myhrf[,2] , type='l' , col='red',pch=2)
points(   myhrf[,3] , type='l' , col='green',pch=2)
dev.off()
