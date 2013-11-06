library( ANTsR )
library( randomForest)
library( vegan )
print("a random forest test on fmri")
subnum<-"009"
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
mat<-( as.matrix( residuals( lm( mat  ~ as.matrix( nuis ) ) ) ) )
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
threshinit<-FALSE
img1<-antsImageClone( mask )
img1[ mask > 0 ]<- myhrf[,1] %*% mat 
if ( threshinit ) img1[ img1 < mean(img1) ] <- 0
img2<-antsImageClone( mask )
img2[ mask > 0 ]<- myhrf[,2] %*% mat 
if ( threshinit ) img2[ img2 < mean(img2) ] <- 0
img3<-antsImageClone( mask )
img3[ mask > 0 ]<- myhrf[,3] %*% mat 
if ( threshinit ) img3[ img3 < mean(img3) ] <- 0
img4<-antsImageClone( mask )
img4[ mask > 0 ]<- myglobsigresid %*% mat
if ( threshinit ) img4[ img4 < mean(img4) ] <- 0
initlist<-list( img1,img2,img3,img4 )
ff<-sparseDecom2( inmatrix=list(mat, cbind(myhrf,myglobsigresid)), inmask=list(mask,NA),
                  perms=0, its=45,
                  mycoption=1, sparseness=c(-0.1,-0.35) ,nvecs=6,
                  smooth=1, cthresh=c(10,0),ell1 = 11 , z=-1 )
#                ,initializationList = initlist )
for ( i in 1:length(ff$eig1) ) {
  visimg<-antsImageClone( ff$eig1[[i]] )
  ImageMath(3,visimg,'abs',visimg)
  ImageMath(3,visimg,'Normalize',visimg)
  plotANTsImage( myantsimage=temimg, functional=list(visimg) , slices="2x25x1", ,axis=0, color=c("red") , threshold="0.01x1", outname=paste("temp",i,".jpg",sep='') ) 
}
print( ff$eig2 )
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
mat2<-( as.matrix( residuals( lm( mat2  ~ as.matrix( nuis2 ) ) ) ) )
myglobsig2<-apply( mat2, FUN=mean, MARGIN=1 )
mat2[,1]<-myglobsig2
mysccanimages<-imageListToMatrix( imageList=ff$eig1, mask=mask)
mysccanpreds <-decostand( mat2 , method="standardize" , MARGIN=2 )  %*% t( mysccanimages )
colnames( mysccanpreds )<-colnames( ff$projections )
mydf2<-data.frame(  imgs = mysccanpreds )

mypred<-predict( my.rf1 , newdata = mydf2 )
print( cor.test(  as.numeric(my.rf1$predicted) , myhrf2[,1] ) )
plot(  myhrf2[,1] , type='l' )
points(    as.numeric(my.rf1$predicted) , type='l' , col='green')

mypred<-predict( my.rf2 , newdata = mydf2 )
print( cor.test(  as.numeric(my.rf2$predicted) , myhrf2[,2] ) )
plot(  myhrf2[,2] , type='l' )
points(    as.numeric(my.rf2$predicted) , type='l' , col='green')

mypred<-predict( my.rf3 , newdata = mydf2 )
print( cor.test(  as.numeric(my.rf3$predicted) , myhrf2[,3] ) )
plot(  myhrf2[,3] , type='l' )
points(    as.numeric(my.rf3$predicted) , type='l' , col='green')
