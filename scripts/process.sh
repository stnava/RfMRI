task=task002
subs="sub007 sub001 sub002 sub003 sub004 sub005 sub006	 sub008 sub009 sub010"
subs=sub002
homedir=`pwd`
for x in $subs ; do
  if [[ ! -s ${homedir}/${x}/BOLD/${task}_run001/avg.nii.gz ]] ; then 
    cd ${homedir}/${x}/BOLD/${task}_run001/
    design=${homedir}/${x}/model/model001/onsets/${task}_run002/cond003.txt
    design=${homedir}/sub004/model/model001/onsets/${task}_run001/cond001.txt 
    ${homedir}/scripts/process_bold.R  -d $design -t 2.5 -n run001 -s sub001
    cd ${homedir}/${x}/BOLD/${task}_run002/
    ${homedir}/scripts/process_bold.R  -d $design -t 2.5 -n run001 -s sub001
  fi
#########################################################
  avg1=${homedir}/${x}/BOLD/${task}_run001/avg.nii.gz 
  avg2=${homedir}/${x}/BOLD/${task}_run002/avg.nii.gz
  bold1=${homedir}/${x}/BOLD/${task}_run001/sccan.nii.gz  
  bold2=${homedir}/${x}/BOLD/${task}_run002/sccan.nii.gz 
  # multivariate comparison
  if [[ -s $avg1 ]] &&  [[ -s $avg2 ]] &&  [[ -s $bold1 ]] &&  [[ -s $bold2 ]] ; then 
  mkdir -p ${homedir}/${x}/BOLD/comparem/
  cd ${homedir}/${x}/BOLD/comparem/
  ${homedir}/scripts/ants_compare.sh $avg1 $avg2 $bold1 $bold2 x
  pre=${homedir}/${x}/BOLD/${task}_run001/
  ${homedir}/scripts/activity_cross_validation.R --boldimg ${pre}mat.mha --mask ${pre}mask.nii.gz --hrf ${pre}antsr_hrf.csv -a bold2to1_x_warped.nii.gz -o ./prediction_x.csv 

  ${homedir}/scripts/ants_compare.sh $avg2 $avg1 $bold2 $bold1 y
  pre=${homedir}/${x}/BOLD/${task}_run002/
  ${homedir}/scripts/activity_cross_validation.R --boldimg ${pre}mat.mha --mask ${pre}mask.nii.gz --hrf ${pre}antsr_hrf.csv -a bold2to1_y_warped.nii.gz -o ./prediction_y.csv 

  else
  echo subject $x failed multivar - check data  
  fi
  # univariate comparison
  bold1=${homedir}/${x}/BOLD/${task}_run001/betast.nii.gz 
  bold2=${homedir}/${x}/BOLD/${task}_run002/betast.nii.gz 
  # some ants stuff 
  if [[ -s $avg1 ]] &&  [[ -s $avg2 ]] &&  [[ -s $bold1 ]] &&  [[ -s $bold2 ]] ; then 
  mkdir -p ${homedir}/${x}/BOLD/compareu/
  cd ${homedir}/${x}/BOLD/compareu/
  ${homedir}/scripts/ants_compare.sh $avg1 $avg2 $bold1 $bold2 x
  pre=${homedir}/${x}/BOLD/${task}_run001/
  ${homedir}/scripts/activity_cross_validation.R --boldimg ${pre}mat.mha --mask ${pre}mask.nii.gz --hrf ${pre}antsr_hrf.csv -a bold2to1_x_warped.nii.gz -o ./prediction_x.csv 

  ${homedir}/scripts/ants_compare.sh $avg2 $avg1 $bold2 $bold1 y
  pre=${homedir}/${x}/BOLD/${task}_run002/
  ${homedir}/scripts/activity_cross_validation.R --boldimg ${pre}mat.mha --mask ${pre}mask.nii.gz --hrf ${pre}antsr_hrf.csv -a bold2to1_y_warped.nii.gz -o ./prediction_y.csv 
  else
  echo subject $x failed univar - check data  
  fi
done 
