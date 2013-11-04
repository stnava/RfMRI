#### group-wise reproducibility of gorgolewski finger-tapping & covert verbs generation ####
ct=0
tstats=( 1.5 2 2.5 3 3.5 4 ) 
myspars="  0.01 0.02 0.04 0.08 0.12 0.16 " 
myspars="  0.16 0.12 0.08 0.04 0.02 0.01 " 
for spar in $myspars ; do 
######
 spar2=${tstats[${ct}]}
######
  ./scripts/process_bold_group.R --tr 2.5 --design task002 \
      --run run001 --templatemask ./template/aal.nii.gz \
      --bold group  --output GroupTask002Run001 --statval ${spar}x${spar2}
######
  ./scripts/process_bold_group.R --tr 2.5 --design task002 \
      --run run002 --templatemask ./template/aal.nii.gz \
      --bold group  --output GroupTask002Run002 --statval ${spar}x${spar2}
######
  ./scripts/process_bold_group.R --tr 2.5 --design task003 \
      --run run001 --templatemask ./template/aal.nii.gz \
      --bold group  --output GroupTask003Run001 --statval ${spar}x${spar2}
######
  ./scripts/process_bold_group.R --tr 2.5 --design task003 \
      --run run002 --templatemask ./template/aal.nii.gz \
      --bold group  --output GroupTask003Run002  --statval ${spar}x${spar2}
######evaluation#######
  for x in 2 3  ; do 
    ThresholdImage 3 GroupTask00${x}Run001sccan.nii.gz temp.nii.gz 1.e-9 1
    ThresholdImage 3 GroupTask00${x}Run002sccan.nii.gz temp2.nii.gz 1.e-9 1
    ImageMath 3 Z${x}_${spar}_eval DiceAndMinDistSum temp.nii.gz temp2.nii.gz
    ThresholdImage 3 GroupTask00${x}Run001group_betasT.nii.gz temp.nii.gz 1.e-9 1.e9
    ThresholdImage 3 GroupTask00${x}Run002group_betasT.nii.gz temp2.nii.gz 1.e-9 1.e9
    ImageMath 3 U${x}_${spar2}_eval DiceAndMinDistSum temp.nii.gz temp2.nii.gz
    rm [U,Z]*mds.nii.gz [U,Z]*dice.nii.gz
  done
######
let ct=${ct}+1
######
done


