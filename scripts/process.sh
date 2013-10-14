
subs="sub001 sub002 sub003 sub004 sub005 sub006	sub008 sub007 sub009 sub010"
homedir=`pwd`
for x in $subs ; do
#  if [[ -s ${homedir}/${x}/BOLD/task003_run001/bold.nii.gz ]] ; then 
    cd ${homedir}/${x}/BOLD/task003_run001/
    design=${homedir}/${x}/model/model001/onsets/task003_run002/cond003.txt
    design=${homedir}/sub004/model/model001/onsets/task003_run001/cond001.txt 
    ${homedir}/scripts/process_bold.R  -d $design -t 2.5 -n run001 -s sub001
    cd ${homedir}/${x}/BOLD/task003_run002/
    ${homedir}/scripts/process_bold.R  -d $design -t 2.5 -n run001 -s sub001
#  fi
  cd ${x}/BOLD/
#########################################################
  avg1=${homedir}/${x}/BOLD/task003_run001/avg.nii.gz 
  avg2=${homedir}/${x}/BOLD/task003_run002/avg.nii.gz
  bold1=${homedir}/${x}/BOLD/task003_run001/sccan.nii.gz 
  bold2=${homedir}/${x}/BOLD/task003_run002/sccan.nii.gz 
  # multivariate comparison
  if [[ -s $avg1 ]] &&  [[ -s $avg2 ]] &&  [[ -s $bold1 ]] &&  [[ -s $bold2 ]] ; then 
  mkdir -p ${homedir}/${x}/BOLD/comparem/
  cd ${homedir}/${x}/BOLD/comparem/
  ln -s $avg1 .
  ln -s $bold1 .
  ${homedir}/scripts/ants_compare.sh $avg1 $avg2 $bold1 $bold2 x
  pre=${homedir}/${x}/BOLD/task003_run001/
  ${homedir}/scripts/activity_cross_validation.R --boldimg ${pre}mat.mha --mask ${pre}mask.nii.gz --hrf ${pre}antsr_hrf.csv -a bold2to1_x_warped.nii.gz -o ./prediction.csv 
  else
  echo subject $x failed multivar - check data  
  fi
  # univariate comparison
  bold1=${homedir}/${x}/BOLD/task003_run001/betast.nii.gz 
  bold2=${homedir}/${x}/BOLD/task003_run002/betast.nii.gz 
  # some ants stuff 
  if [[ -s $avg1 ]] &&  [[ -s $avg2 ]] &&  [[ -s $bold1 ]] &&  [[ -s $bold2 ]] ; then 
  mkdir -p ${homedir}/${x}/BOLD/compareu/
  cd ${homedir}/${x}/BOLD/compareu/
  ln -s $avg1 .
  ln -s $bold1 .
  ${homedir}/scripts/ants_compare.sh $avg1 $avg2 $bold1 $bold2 x
  pre=${homedir}/${x}/BOLD/task003_run001/
  ${homedir}/scripts/activity_cross_validation.R --boldimg ${pre}mat.mha --mask ${pre}mask.nii.gz --hrf ${pre}antsr_hrf.csv -a bold2to1_x_warped.nii.gz -o ./prediction.csv 
  else
  echo subject $x failed univar - check data  
  fi
done 
