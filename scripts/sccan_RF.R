library( ANTsR )
library( randomForest)
library( vegan )
takeoutresid <- TRUE
print(paste("a random forest test on fmri ... takeoutresid? ",takeoutresid))
subnum<-"010"
setwd(paste("/Users/stnava/data/data_gorgolewski/RfMRI/group_analysis/sub",subnum,"/task003/run001",sep=''))
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
mat<-timeseries2matrix( fmri , mask )
myglobsig<-apply( mat, FUN=mean, MARGIN=1 )
boldthresh<-1200 # mean(myglobsig)-6*sd(myglobsig)
selector<-myglobsig > boldthresh
selector[1:(nbadframes-1)] <- FALSE
mat<-subset( mat , selector )
nuis<-subset( nuis , selector )
myhrf<-subset( ohrf , selector )
if ( takeoutresid )  mat<-( as.matrix( residuals( lm( mat  ~ as.matrix( nuis ) ) ) ) )
myglobsig<-apply( mat, FUN=mean, MARGIN=1 )
mat[,1]<-myglobsig
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
threshinit<-3
img1<-antsImageClone( mask )
img1[ mask > 0 ]<- myhrf[,1] %*% mat 
if ( threshinit > 0 ) img1[ img1 <  (mean(img1)+threshinit*sd(img1)) ] <- 0
img2<-antsImageClone( mask )
img2[ mask > 0 ]<- myhrf[,2] %*% mat 
if ( threshinit > 0 ) img2[ img2 <  (mean(img2)+threshinit*sd(img2)) ] <- 0
img3<-antsImageClone( mask )
img3[ mask > 0 ]<- myhrf[,3] %*% mat 
if ( threshinit > 0 ) img3[ img3 <  (mean(img3)+threshinit*sd(img3)) ] <- 0
img4<-antsImageClone( mask )
img4[ mask > 0 ]<- myglobsigresid %*% mat
if ( threshinit > 0) img4[ img4 < (mean(img4)+threshinit*sd(img4))  ] <- 0
img5<-antsImageClone( mask )
img5[ mask > 0 ]<- svd( mat )$u[5,] %*% mat
if ( threshinit > 0 ) img5[ img5 < (mean(img5)+threshinit*sd(img5)) ] <- 0
initlist<-list( img1,img2,img3,img4, img5 )
ff<-sparseDecom2( inmatrix=list(mat, cbind(myhrf,myglobsigresid)), inmask=list(mask,NA), perms=0, its=22,mycoption=1, sparseness=c( -0.1, -0.2 ) , nvecs= length( initlist ) , smooth=1, cthresh=c(10,0), ell1 = 11 , z=-1 , initializationList=initlist )
print( ff$eig2 )
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
mydf<-data.frame( binpreds = myhrf[,1] , imgs = ff$projections )
my.rf1<-randomForest( binpreds ~ . , data=mydf )
mydf<-data.frame( binpreds = myhrf[,2] , imgs = ff$projections )
my.rf2<-randomForest( binpreds ~ . , data=mydf )
mydf<-data.frame( binpreds = myhrf[,3] , imgs = ff$projections )
my.rf3<-randomForest( binpreds ~ . , data=mydf )

######################################################################
fmri2<-antsImageRead( paste('../run002/sub',subnum,'_group.nii.gz',sep='') , 4 )
nuis2<-read.csv('nuis.csv')
mat2<-timeseries2matrix( fmri2 , mask )
myglobsig2<-apply( mat2, FUN=mean, MARGIN=1 )
selector2<-myglobsig2 > boldthresh
selector2[1:(nbadframes-1)] <- FALSE
mat2<-subset( mat2 , selector )
nuis2<-subset( nuis2 , selector )
myhrf2<-subset( ohrf , selector )
if ( takeoutresid ) mat2<-( as.matrix( residuals( lm( mat2  ~ as.matrix( nuis2 ) ) ) ) )
myglobsig2<-apply( mat2, FUN=mean, MARGIN=1 )
mat2[,1]<-myglobsig2
mysccanimages<-imageListToMatrix( imageList=ff$eig1, mask=mask)
mysccanpreds <-decostand( mat2 , method="standardize" , MARGIN=2 )  %*% t( mysccanimages )
colnames( mysccanpreds )<-colnames( ff$projections )
mydf2<-data.frame(  imgs = mysccanpreds )
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
