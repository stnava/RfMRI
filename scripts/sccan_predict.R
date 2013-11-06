library(ANTsR)
library( randomForest )
library( vegan )
off<-( 2 ) 
mysub<-"sub003"
myrun<-"run001"
myrun2<-"run002"
tasknum<-2
if ( tasknum == 3 ) domotor<-TRUE else domotor <- FALSE 
mytask<-paste("task00",tasknum,sep='')
mybasefn<-paste(".//group_analysis",mysub,mytask,myrun,"",sep='/')
mybasefn2<-paste(".//group_analysis",mysub,mytask,myrun2,"",sep='/')
fn1<-paste(mybasefn,mysub,"_group.nii.gz",sep='')
fn2<-paste(mybasefn2,mysub,"_group.nii.gz",sep='')
mask<-antsImageRead( "./template/aal.nii.gz" , 3 )
mask2<-antsImageRead(paste("OptGroupTask00",tasknum,"Run001sccan.nii.gz",sep=''),3)
th<-c(1,90)
# mask[ mask > th[2] | mask < th[1]  | mask2 < 1.e-6  ]<-0
# mask[ mask >= th[1] & mask <= th[2] & mask2 > 1.e-6 ]<-1
mask[ mask > th[2] | mask < th[1]   ]<-0
mask[ mask >= th[1] & mask <= th[2] ]<-1
mask2vec<-mask2[ mask > 0 ]
mask2vec<-mask2vec/sum(mask2vec)
nu1<-read.csv(paste(mybasefn,"nuis.csv",sep=''))[,c(2:9)]
nu2<-read.csv(paste(mybasefn2,"nuis.csv",sep=''))[,c(2:9)]
fmri1<-timeseries2matrix( antsImageRead( fn1, 4 ) , mask ) 
fmri2<-timeseries2matrix( antsImageRead( fn2, 4 ) , mask ) 
# fmri1<-decostand( timeseries2matrix( antsImageRead( fn1, 4 ) , mask ) , method="standardize" , MARGIN=2 )
# fmri2<-decostand( timeseries2matrix( antsImageRead( fn2, 4 ) , mask ) , method="standardize" , MARGIN=2 )
hrf1<-read.csv(paste(mybasefn,"antsr_hrf.csv",sep=''))
##############################################################
mytemin<-which( hrf1[,1] > 0 )
mytem<-mytemin[ mytemin > 25 & mytemin < 44 ]
if ( domotor ) mytem<-mytemin[ mytemin > 40 & mytemin < 77 ]
mytem<-c( (mytem[1]-off):(mytem[1]-1),mytem,c( (mytem[length(mytem)]+1):(mytem[length(mytem)]+off) ) )
mytem2<-mytemin[ mytemin > 100 & mytemin < 120 ]
if ( domotor ) mytem2<-mytemin[ mytemin > 111 & mytemin < 128 ]
mytem2<-c( (mytem2[1]-off):(mytem2[1]-1),mytem2,c( (mytem2[length(mytem2)]+1):(mytem2[length(mytem2)]+off) ) )
##############################################################
rfmri1<-( as.matrix( residuals( lm( fmri1  ~ as.matrix( nu1 ) ) ) ) )
rfmri2<-( as.matrix( residuals( lm( fmri2  ~ as.matrix( nu2 ) ) ) ) )
template<-whiten( rfmri1[ mytem , ] ) + whiten( rfmri1[ mytem2 , ] )
# dev.new()
ntests<-( nrow(fmri1)-nrow(template) )
myccas<-rep(NA,ntests)
ncomp<-2
refsvd<-svd( template )$u[,2:ncomp]
templateproj<- template %*% mask2vec 
for ( x in (off+1):ntests )
  {
  subfmri<-rfmri2[c((x-off):(x+nrow(template)-off-1)), ]
#  locsvd<-svd( subfmri )$u[,2:ncomp]
  myccas[x]<-( cor.test( templateproj, subfmri %*% mask2vec )$est  )*( 1 ) 
#  sccan<-sparseDecom2( inmatrix=list( template, subfmri ) , nvecs=2, robust=0,inmask = c( mask, mask),     its=3, mycoption=1 ,  perms=1, sparseness=c( -1 , -1 ), z=-1, ell1=11 , smooth=0, cthresh=c(0,0))
#  myccas[x]<-sccan$ccasummary[2,2]
  print( paste( x/ntests*100.,"%" ))
  yv<-0.0
  plot( myccas, type='l', ylim=c(yv,1))
  points( ( hrf1[,1] - shift(hrf1[,1],1) ) > 0.9 ,type='l',col='red')
  }
