#!/usr/bin/env Rscript
smth<-1
nrandrun<-4
ndvis<-8
zval<-( -1 )
clustval<-25
timefilter<-FALSE
residdesign<-FALSE
CORAC<-FALSE
dounivar<-TRUE 
sparval <- 0.2
sparval2<- -0.3
domotor <- FALSE
library(getopt)
initlist<-list()
options(digits=3)
Args <- commandArgs()
self<-Args[4]
self<-substring(self,8,nchar(as.character(self)))
spec = c(
  'templatemask', 's', 1, "character" ," 3D template mask name ",
  'design', 'd', 1, "character" ," design matrix ",
  'bold', 'b', 1, "character" ," BOLD image ",
  'tr', 't', 1, "character" ," BOLD tr ",
  'smoother', 'm', 0, "character" ," BOLD spatial smoothing ",
  'run', 'n', 1, "character" ," which run ",
  'help'     , 'h', 0, "logical" ," print the help ", 
  'statval', 'a', 1, "character" ," sparseness for multivariate or t-tstat for univariate - of form 0.1x3.0 where the 2nd is for unviarate",
  'output', 'o', 1, "character" ," output ")
spec=matrix(spec,ncol=5,byrow=TRUE)
opt = getopt(spec)
if ( !is.null(opt$help) || length(opt) == 1 ) {
#print a friendly message and exit with a non-zero error code
cat("\n")
cat(paste(self,"\n"))
for ( x in 1:nrow(spec) ) {
cat("\n")
longopt<-paste("--",spec[x,1],sep='')
shortopt<-paste("-",spec[x,2],sep='')
hlist<-paste(shortopt,"|",longopt,spec[x,5],"\n \n")
# print(hlist,quote=F)
cat(format(hlist, width=40, justify = c("left")))
}
q(status=1);
}

#../scripts/process_bold_group.R --tr 2.5 --design task002 --run run002
#  --templatemask ./template/aal.nii.gz  --bold group  --output W ;  
if ( is.null( opt$tr ) ) {
  print( 'no tr , automating parameters or use --tr option')
  opt$tr<-2.5
  opt$statval<-as.character( "0.1x3.5" )
  opt$design<-"task002"
  opt$run<-"run001"
  opt$templatemask<-"./template/aal.nii.gz"
  opt$bold<-"group" 
  opt$output<-"TEST"
}
if ( is.null( opt$statval ) ) {
  opt$statval<-"0.1x3.5" 
}
mysplit<-strsplit(opt$statval,as.character("x"))
sparval<-as.numeric( mysplit[[1]][1] )
statval<-as.numeric( mysplit[[1]][2] )
if ( is.null( opt$design ) ) {
  print( 'no design , quitting , use --design option')
  q()
}
if ( is.null( opt$run ) ) {
  print( 'no run , quitting , use --run option')
  q()
}
library(ANTsR)
library(lme4)
library(nlme)
if ( is.null( opt$templatemask ) ) {
  print( 'no templatemask , quitting , use --templatemask option')
  q()
} else {
  maskin<-antsImageRead( opt$templatemask , 3 )
  mask<-antsImageClone( maskin )
  th<-c(1,90)
  mask[ mask > th[2] | mask < th[1] ]<-0
  mask[ mask >= th[1] & mask <= th[2] ]<-1
  print( paste( "mask size", sum( mask > 0  ) ))
}
#######################################################################
if ( is.null( opt$bold ) ) {
  print( 'no bold image filename, quitting, use --bold option')
  q()
} else {
  searchstring<-paste("group_analysis/sub*/",opt$design,"/",opt$run,"/*group.nii.gz",sep='')
  fns<-Sys.glob(  searchstring )
  searchstring<-paste("group_analysis/sub*/",opt$design,"/",opt$run,"/nuis.csv",sep='')
  fnn<-Sys.glob(  searchstring )
  if ( length( fns ) != length( fnn ) |  length( fnn )  == 1  )
    {
    print( paste("You have ",length(fns)," subjects ... fewer than 3 not recommended.  quitting if only 1.") )
    print( paste("N-nuisance files",length( fnn )," N-bold_group.nii.gz files ", length( fns ),"." ) )
    print( paste("OR -- Number of nuisance files does not match number of *group.nii.gz files => quitting.") )
    q()
    }
  print( fns )
  print( fnn )
  fmri<-antsImageRead( fns[1]  ,4)
  ImageMath(4,fmri,'SliceTimingCorrection',fmri,0)
  smoother<-0
  if ( ! is.null( opt$smoother ) ) {
    smoother <- as.numeric( opt$smoother )
  }
  if ( smoother > 0.001 ) SmoothImage(4,fmri,smoother,fmri)
  dowhite<-T
  mat<-timeseries2matrix( fmri, mask )
  if ( opt$design == "task003" ) domotor<-TRUE
  print(paste("Do Motor",domotor,opt$design," run ", opt$run ))
  tr<-as.numeric( opt$tr )
  nbadframes<-5
  blockfing = c(0, 36, 72, 108,144)+nbadframes
  blockfoot <- blockfing + 12
  blockmout <- blockfoot + 12 
  blocko = c(1, 24, 48, 72, 96, 120, 144 )+nbadframes # covert 
  if ( domotor ) ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfing , durations=rep(  6,  length( blockfing ) ) ,  rt=tr ) else ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blocko , durations=rep(  12,  length( blocko ) ) ,  rt=tr )
  if ( domotor ) {
    ohrf2 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfoot , durations=rep(  6,  length( blockfoot ) ) ,  rt=tr )
    ohrf3 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockmout , durations=rep(  6,  length( blockmout ) ) ,  rt=tr )
  }
  if ( domotor ) ohrf<-cbind( ohrf , ohrf2, ohrf3 ) else ohrf<-cbind( ohrf )
  ohrf[1:4,]<-0 # first few frames are junk
  nuis<-read.csv(fnn[1])
  bnuis<-as.matrix( data.frame( globalsignal = nuis$myvarsin.globalsignal, motion1=nuis$motion1, motion2=nuis$motion2, motion3=nuis$motion3, compcorr1=nuis$compcorr1, compcorr2=nuis$compcorr2, compcorr3=nuis$compcorr3  ) )
  bnuis<-bnuis[nbadframes:nrow(mat), ]
  ohrf<-ohrf[nbadframes:nrow(mat), ]
  mat<-mat[nbadframes:nrow(mat), ]
  myglobsig<-apply( mat, FUN=mean, MARGIN=1 )
  boldthresh<-1200 # mean(myglobsig)-6*sd(myglobsig)
  mat<-subset( mat , myglobsig > boldthresh )
  bnuis<-subset( bnuis , myglobsig > boldthresh )
  ohrf<-as.matrix( subset( ohrf , myglobsig > boldthresh ) )
  myglobsig<-apply( mat, FUN=mean, MARGIN=1 )
  plot( myglobsig , type='l' )
  if ( dowhite ) mat<-whiten( mat ) 
  bhrf<-ohrf
  print("OBHRF")
  print( dim( bhrf ) )
  subjid<-rep(1,nrow(mat) )
  whichsubs<-2:length(fns)
  for ( i in whichsubs )
    {
    print( fns[i] )
    fmri<-antsImageRead( fns[i]  ,4)
    ImageMath(4,fmri,'SliceTimingCorrection',fmri,0)
    if ( smoother > 0.001 ) SmoothImage(4,fmri,smoother,fmri)
    locmat<-timeseries2matrix( fmri, mask )
    if ( domotor ) ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfing , durations=rep(  12,  length( blockfing ) ) ,  rt=tr ) else ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blocko , durations=rep(  12,  length( blocko ) ) ,  rt=tr ) 
    if ( domotor ) {
      ohrf2 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfoot , durations=rep(  12,  length( blockfoot ) ) ,  rt=tr )
      ohrf3 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockmout , durations=rep(  12,  length( blockmout ) ) ,  rt=tr )
    }
    if ( dmotor )  ohrf<-cbind( ohrf , ohrf2, ohrf3 ) else ohrf<-cbind( ohrf )
    print( paste("OHRF",domotor ) )
    print( dim(ohrf) )
    nuis<-read.csv(fnn[i]) 
    nuis<-as.matrix( data.frame( globalsignal = nuis$myvarsin.globalsignal, motion1=nuis$motion1, motion2=nuis$motion2, motion3=nuis$motion3, compcorr1=nuis$compcorr1, compcorr2=nuis$compcorr2, compcorr3=nuis$compcorr3  ) )
    ohrf<-ohrf[nbadframes:nrow(locmat),]
    nuis<-nuis[nbadframes:nrow(locmat),]
    locmat<-locmat[nbadframes:nrow(locmat),]
    myglobsig<-apply( locmat, FUN=mean, MARGIN=1 )
    locmat<-subset( locmat , myglobsig > boldthresh )
    nuis<-subset( nuis , myglobsig > boldthresh )
    ohrf<-as.matrix( subset( ohrf , myglobsig > boldthresh ) )
    print( paste( "LOCMAT DIM" ,dim(locmat) ) )
    myglobsig2<-apply( locmat, FUN=mean, MARGIN=1 )
    plot( myglobsig2 , type='l' )
    if ( dowhite ) locmat<-whiten( locmat )  
    subjid<-c( subjid, rep(i,nrow(locmat)) )
    mat<-rbind( mat , locmat )
    print( paste( "BHRF DIM" ,dim(bhrf) ) )
    print( paste( "OHRF DIM" ,dim(ohrf) ) )
    bhrf<-rbind( bhrf, ohrf )
    bnuis<-rbind( bnuis , nuis )
  }
} 
betas<-rep(NA, ncol(mat) )
pvals<-betas
bnuis<-data.frame( bnuis )
print( names( bnuis ) )
if ( TRUE ) {
progress <- txtProgressBar( min = 0, max = length(betas), style = 3 )
for ( i in 1:length(betas) )
  {
  vox<-mat[ , i ] 
  arval<-ar(vox,FALSE,2)$ar
  arval<-shift(vox,1)*arval[1]+shift(vox,2)*arval[2]
  mdl<-lm( vox ~  bhrf  + motion1 + motion2 + motion3 + compcorr1 + compcorr2 + compcorr3  + globalsignal + subjid , data = bnuis )
  betas[i]<-coefficients( summary( mdl ) )[2,3]
#  mdl <- lme(vox ~  bhrf + motion1 + motion2 + motion3 + compcorr1 + compcorr2 + compcorr3 + globalsignal , random = ( ~ 1 | subjid) , data = bnuis )
#  betas[i]<-summary(mdl)$tTable[2,4]
  pvals[i]<-2*pt(-abs(betas[i]),df= (nrow(mat)-1) )
  if ( i %% 5000 == 0 | i == 1 )
    {
    print(summary(mdl))
    print( paste( betas[i] , max( abs( betas ) , na.rm = T ), min( abs( pvals ) , na.rm = T ) ) )
    }
  setTxtProgressBar(progress, i)
  }
close(progress)
print("done") 
print( max( betas ) )
print( min( betas ) )
print(paste("qvals:", min( p.adjust( pvals ) ) ) )
betaimg<-antsImageClone( mask )
betaimg[ mask > 0 ]<-betas
antsImageWrite(betaimg,paste(opt$output,"group_betas.nii.gz",sep=''))
betaimgt<-antsImageClone( mask )
betaimgt[ mask > 0 ]<-0
betaimgt[ betaimg > statval ]<-betaimg[ betaimg > statval ]
antsImageWrite(betaimgt,paste(opt$output,"group_betasT.nii.gz",sep=''))
}
################################################
##################SPARSE-CCA####################
################################################
print("BeginResid")
rmat<-as.matrix( residuals( lm( as.matrix(mat) ~ 1 + motion1 + motion2 + motion3 + compcorr1 + compcorr2 + compcorr3 + globalsignal + as.factor(subjid) , data = bnuis ) ) )
print("EndResid")
mypreds<-as.matrix( cbind( bhrf , as.numeric(  bhrf[,1] ) ) )
nv<-ncol(mypreds)
print(paste("Get",nv,"variates"))
sccan<-sparseDecom2( inmatrix=list( rmat , mypreds ), inmask = c( mask , NA ) ,
                     sparseness=c( sparval, 0.1 ), nvecs=nv, its=25, smooth=1,
                     perms=0, cthresh = c(10, 0) , robust=0, mycoption=1,
                     z=-1    , ell1=11 )
print( sccan$eig2 )
for ( ee in 1:length( sccan$eig1 ) ) ImageMath(3,sccan$eig1[[ee]],"abs",sccan$eig1[[ee]]) 
eigres<-sccan$eig1[[1]]
cor1<-abs( cor.test( bhrf[,1] , sccan$projections[,1] )$est )
cor2<-abs( cor.test( bhrf[,1] , sccan$projections[,2] )$est )
print( paste("cor1",cor1  ))
for ( ss in unique( subjid ) )
  {
    corsN<-abs( cor.test( bhrf[ subjid == ss ,1] , sccan$projections[ subjid == ss ,1] )$est )
    print( paste(ss,"corsN",corsN  ))
  }
print( paste("cor2",cor2  ))
if ( cor2 > cor1 ) {
  eigres<-sccan$eig1[[2]]
  for ( ss in unique( subjid ) )
    {
    corsN<-abs( cor.test( bhrf[ subjid == ss ,1] , sccan$projections[ subjid == ss ,2] )$est )
    print( paste(ss,"corsN",corsN  ))
    }
}
if ( length(sccan$eig1) > 2 )
  {
  cor3<-abs( cor.test( bhrf[,1] , sccan$projections[,3] )$est )
  print( paste("cor3",cor3  ))
  if ( cor3 > max( c( cor1, cor2 ) ) ) {
    eigres<-sccan$eig1[[3]]
    for ( ss in unique( subjid ) )
      {
        corsN<-abs( cor.test( bhrf[ subjid == ss ,1] , sccan$projections[ subjid == ss ,3] )$est )
        print( paste(ss,"corsN",corsN  ))
      }
    }
  }
if ( length(sccan$eig1) > 3 )
  {
  cor4<-abs(cor.test( bhrf[,1] , sccan$projections[,4] )$est )
  print( paste("cor4",cor4  ))
  if ( cor4 > max( c( cor1, cor2, cor3 ) ) ) eigres<-sccan$eig1[[4]]
  }
antsImageWrite( eigres ,  paste(opt$output,"sccan.nii.gz",sep="")  )
antsImageWrite( as.antsImage( rmat ) , paste(opt$output,"rmat.mha",sep="") )
write.csv( bhrf, paste(opt$output,"bhrf.csv",sep="") , quote=F , row.names=F )
