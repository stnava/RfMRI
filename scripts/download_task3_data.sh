
for x in   sub003   ; do 

  mkdir ${x}
  if [[ ! -s ./${x}/BOLD/task003_run001/bold.nii.gz ]] ; then 
    wget -r ftp://anonymous:password@climb.genomics.cn:/pub/10.5524/100001_101000/100051/data/data/${x}/BOLD/task003_run001/bold.nii.gz
  fi 
  if [[ ! -s ./${x}/BOLD/task003_run001/bold.nii.gz ]] ; then 
    wget -r ftp://anonymous:password@climb.genomics.cn:/pub/10.5524/100001_101000/100051/data/data/${x}/BOLD/task003_run002/bold.nii.gz
  fi 
  mv .//climb.genomics.cn/pub/10.5524/100001_101000/100051/data/data/${x}/BOLD ./${x}/ 

done
`
