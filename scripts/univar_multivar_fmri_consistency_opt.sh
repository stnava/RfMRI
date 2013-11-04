#### group-wise reproducibility of gorgolewski finger-tapping & covert verbs generation ####
#### with optimal parameters for each strategy / task ####
######
  moo=" --smoother 1.5  "
  ./scripts/process_bold_group.R --tr 2.5 --design task002 \
      --run run001 --templatemask ./template/aal.nii.gz \
      --bold group  --output OptGroupTask002Run001 --statval 0.12x2.0 $moo
######
  ./scripts/process_bold_group.R --tr 2.5 --design task002 \
      --run run002 --templatemask ./template/aal.nii.gz \
      --bold group  --output OptGroupTask002Run002 --statval 0.12x2.0 $moo
######
  ./scripts/process_bold_group.R --tr 2.5 --design task003 \
      --run run001 --templatemask ./template/aal.nii.gz \
      --bold group  --output OptGroupTask003Run001 --statval 0.02x3.5 $moo
######
  ./scripts/process_bold_group.R --tr 2.5 --design task003 \
      --run run002 --templatemask ./template/aal.nii.gz \
      --bold group  --output OptGroupTask003Run002  --statval 0.02x3.5 $moo 
######evaluation#######
  for x in 2 3  ; do 
    ThresholdImage 3 OptGroupTask00${x}Run001sccan.nii.gz temp.nii.gz 1.e-9 1
    ThresholdImage 3 OptGroupTask00${x}Run002sccan.nii.gz temp2.nii.gz 1.e-9 1
    ImageMath 3 OptZ${x}_eval DiceAndMinDistSum temp.nii.gz temp2.nii.gz
    ThresholdImage 3 OptGroupTask00${x}Run001group_betasT.nii.gz temp.nii.gz 1.e-9 1.e9
    ThresholdImage 3 OptGroupTask00${x}Run002group_betasT.nii.gz temp2.nii.gz 1.e-9 1.e9
    ImageMath 3 OptU${x}_eval DiceAndMinDistSum temp.nii.gz temp2.nii.gz
    rm [OptU,OptZ]*mds.nii.gz [OptU,OptZ]*dice.nii.gz
  done

