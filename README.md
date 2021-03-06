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
# this script uses wget which you install via homebrew on osx
# on linux wget is native or its alternative, yum 
homedir=${PWD}  # root RfMRI directory 
${homedir}/scripts/download_task_data.sh
# now process each subject to group space 
tem=${homedir}/template/template.nii.gz
for task in task002 task003  ; do 
  for myrun in run001 run002 ; do 
    for x in sub001 sub002 sub003 sub004 sub005 sub006 sub007  sub008 sub009 sub010  ; do 
      cd $homedir 
      mysub=${homedir}/${x}/BOLD/${task}_${myrun}/avg.nii.gz
      mybold=${homedir}/${x}/BOLD/${task}_${myrun}/bold.nii.gz
      if [[ -s $mybold ]] ; then 
        antsMotionCorr -d 3 -a $mybold -o $mysub # average bold 
	gdir=${homedir}/group_analysis/${x}/${task}/${myrun}
        if [[ ! -s ${gdir}/${x}_group.nii.gz ]] ; then 
          mkdir -p $gdir 
          ${homedir}/scripts/ants_2_template.sh $tem $mysub $mybold ${gdir}/${x}_group  
          cd ${gdir}
          ${homedir}/scripts/process_bold.R  -d $task -t 2.5 -n ${myrun} -s $x --bold ${x}_group.nii.gz
        fi
      fi 
    done					 
  done
done
```
Now we have first-level statistics done and a map to template space.  

Next do group level statistics.

```
# then for each task / run 
./scripts/process_bold_group.R --tr 2.5 --design task002 \
  --run run001 --templatemask ./template/aal.nii.gz \
  --bold group  --output GroupTask002Run001
```

see the pdf presentation and scripts directory for processing description.

several key-points:  groupwise analysis, nuisance variables, look at globalsignal for subset-selection, random effects, slice-timing, hrf-model selection, whitening ...
