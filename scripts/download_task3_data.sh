task=task002
for x in   sub003   ; do 
  mkdir ${x}
  if [[ ! -s ./${x}/BOLD/${task}_run001/bold.nii.gz ]] ; then 
    wget -r ftp://anonymous:password@climb.genomics.cn:/pub/10.5524/100001_101000/100051/data/data/${x}/BOLD/${task}_run001/bold.nii.gz
  fi 
  if [[ ! -s ./${x}/BOLD/${task}_run002/bold.nii.gz ]] ; then 
    wget -r ftp://anonymous:password@climb.genomics.cn:/pub/10.5524/100001_101000/100051/data/data/${x}/BOLD/${task}_run002/bold.nii.gz
  fi 
  mkdir -p ./${x}/BOLD/${task}_run001/  ./${x}/BOLD/${task}_run002/
  mv .//climb.genomics.cn/pub/10.5524/100001_101000/100051/data/data/${x}/BOLD/${task}_run001/bold.nii.gz ./${x}/BOLD/${task}_run001/bold.nii.gz 
  mv .//climb.genomics.cn/pub/10.5524/100001_101000/100051/data/data/${x}/BOLD/${task}_run002/bold.nii.gz ./${x}/BOLD/${task}_run002/bold.nii.gz 
done
