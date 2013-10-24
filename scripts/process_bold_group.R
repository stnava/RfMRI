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


if ( is.null( opt$bold ) ) {
  print( 'no bold image filename, quitting, use --bold option')
  q()
} else {
  fns<-Sys.glob( opt$bold ) # "*/*/*/*group.nii.gz"
  fns<-Sys.glob( "*/*/*/*group.nii.gz" )
  fmri<-antsImageRead( fns[1]  ,4)
  smoother<-3
  SmoothImage(4,fmri,smoother,fmri)
  mat<-timeseries2matrix( fmri, mask )
  print( dim(mat) )
if ( opt$design == "task003" ) domotor<-TRUE
print(paste("Do Motor",domotor,opt$design))
tr<-as.numeric( opt$tr )
# onsets<-round(stim/tr)
blockfing = c(0, 36, 72, 108,144)
blockfoot <- blockfing + 12
blockmout <- blockfoot + 12 ; blockmout<-blockmout[1:(length(blockmout)-1)]
blocko = c(1, 24, 48, 72, 96, 120, 144 ) # covert 
if ( domotor ) ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfing , durations=rep(  12,  length( blockfing ) ) ,  rt=tr ) else ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blocko , durations=rep(  12,  length( blocko ) ) ,  rt=tr ) 
ohrf2 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfoot , durations=rep(  12,  length( blockfoot ) ) ,  rt=tr )
ohrf3 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockmout , durations=rep(  12,  length( blockmout ) ) ,  rt=tr )
if ( domotor ) ohrf<-cbind( ohrf , ohrf2, ohrf3 ) else ohrf<-cbind( ohrf )
ohrf[1:4,]<-0 # first few frames are junk
bhrf<-ohrf
subjid<-rep(1,nrow(mat) )
  for ( i in 2:length(fns) )
    {
    fmri<-antsImageRead( fns[i]  ,4)
    SmoothImage(4,fmri,smoother,fmri)
    locmat<-timeseries2matrix( fmri, mask )
    subjid<-c( subjid, rep(i,nrow(locmat)) )
    mat<-rbind( mat , locmat )
    if ( domotor ) ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfing , durations=rep(  12,  length( blockfing ) ) ,  rt=tr ) else ohrf <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blocko , durations=rep(  12,  length( blocko ) ) ,  rt=tr ) 
    ohrf2 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockfoot , durations=rep(  12,  length( blockfoot ) ) ,  rt=tr )
    ohrf3 <- hemodynamicRF( scans=dim(fmri)[4] , onsets=blockmout , durations=rep(  12,  length( blockmout ) ) ,  rt=tr )
    if ( domotor ) ohrf<-cbind( ohrf , ohrf2, ohrf3 ) else ohrf<-cbind( ohrf )
    ohrf[1:4,]<-0 # first few frames are junk
    bhrf<-rbind( bhrf, ohrf )
    }
  fns<-Sys.glob( "*/*/*/nuis.csv" )
  nuis<-read.csv(fns[1])
  bnuis<-as.matrix( data.frame( globalsignal = nuis$myvarsin.globalsignal, motion1=nuis$motion1, motion2=nuis$motion2, motion3=nuis$motion3, compcorr1=nuis$compcorr1, compcorr2=nuis$compcorr2, compcorr3=nuis$compcorr3  ) )
  for ( i in 2:length(fns) )
    {
    nuis<-read.csv(fns[i]) 
    nuis<-as.matrix( data.frame( globalsignal = nuis$myvarsin.globalsignal, motion1=nuis$motion1, motion2=nuis$motion2, motion3=nuis$motion3, compcorr1=nuis$compcorr1, compcorr2=nuis$compcorr2, compcorr3=nuis$compcorr3  ) )
    bnuis<-rbind( bnuis , nuis )
    }
} 
betas<-rep(NA, ncol(mat) )
bnuis<-data.frame( bnuis )
print( names( bnuis ) )
progress <- txtProgressBar( min = 0, max = length(betas), style = 3 )
for ( i in 1:length(betas) )
  {
  vox<-mat[ , i ] # + globalsignal
#  mdl<-lm( vox ~  bhrf + motion1 + motion2 + motion3 + compcorr1 + compcorr2 + compcorr3  + subjid , data = bnuis )
  mdl <- lme(vox ~  bhrf + motion1 + motion2 + motion3 + compcorr1 + compcorr2 + compcorr3 + globalsignal , random = ( ~ 1 | subjid) , data = bnuis )
  betas[i]<-summary(mdl)$tTable[2,4]
  if ( i == 1 )
    {
    print(summary(mdl))
    print( betas[i] )
    }
  setTxtProgressBar(progress, i)
  }
close(progress)
print("done") 
print( max( betas ) )
print( min( betas ) )
print( length( betas ) )
print( sum( mask > 0 ) )
mask[ mask > 0 ]<-betas
antsImageWrite(mask,"group_betas.nii.gz")

