RfMRI
=====

antsr reproducibility study of fmri

[slides are here](https://github.com/stnava/RfMRI/raw/gh-pages/rfmri.pdf?raw)

analyze data from:

Single subject fMRI test–retest reliability metrics and confounding factors,
Krzysztof J. Gorgolewski, Amos J. Storkey c, Mark E. Bastin b, Ian Whittle d, Cyril Pernet

[article](http://www.gigasciencejournal.com/content/2/1/6)

by running an organization script then a processing script 

```
./scripts/download_task_data.sh
./scripts/process_bold_group.R
```

see the pdf presentation and scripts directory for processing description.

several key-points:  groupwise analysis, nuisance variables, look at globalsignal for subset-selection, random effects, slice-timing, hrf-model selection, whitening ...
