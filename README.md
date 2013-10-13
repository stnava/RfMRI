RfMRI
=====

antsr reproducibility study of fmri

get data from:


by 

for x in  sub003 sub004 sub005 sub006 sub007  ; do 
  mkdir ${x}
  wget -r ftp://anonymous:password@climb.genomics.cn:/pub/10.5524/100001_101000/100051/data/data/${x}/*
done


then see the scripts directory for processing
