
ct=0
tstats=( 2 3 4 5 6 8 ) 
for spar in 0.01 0.02 0.04 0.08 0.12 0.16 ; do 
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
  ImageMath 3 Z${x}_${spar}_${spar2}_eval DiceAndMinDistSum temp.nii.gz temp2.nii.gz
done
######
 let ct=${ct}+1
######
 done


