RfMRI
=====

antsr reproducibility study of fmri

[slides are here](https://github.com/stnava/RfMRI/raw/gh-pages/rfmri.pdf?raw)

get data from:

Single subject fMRI test–retest reliability metrics and confounding factors,
Krzysztof J. Gorgolewski, Amos J. Storkey c, Mark E. Bastin b, Ian Whittle d, Cyril Pernet

[article](http://www.gigasciencejournal.com/content/2/1/6)

by 

```
for x in  sub003 sub004 sub005 sub006 sub007  ; do 

  mkdir ${x}

  wget -r ftp://anonymous:password@climb.genomics.cn:/pub/10.5524/100001_101000/100051/data/data/${x}/*

  mv .//climb.genomics.cn/pub/10.5524/100001_101000/100051/data/data/${x} . 

done
```



then see the pdf presentation and scripts directory for processing description


several key-points:  groupwise analysis, nuisance variables, look at globalsignal for subset-selection, random effects, slice-timing, hrf-model selection, whitening ...

