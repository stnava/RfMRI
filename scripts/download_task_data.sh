for task in task002 task003 ; do 
  for x in  sub001 sub002 sub003 sub004 sub005 sub006 sub007  sub008 sub009 sub010   ; do 
    for myrun in run001 run002 ; do 
      if [[ ! -s ./${x}/BOLD/${task}_${myrun}/bold.nii.gz ]] ; then 
        wget -r ftp://anonymous:password@climb.genomics.cn:/pub/10.5524/100001_101000/100051/data/data/${x}/BOLD/${task}_${myrun}/bold.nii.gz
	mkdir -p ./${x}/BOLD/${task}_${myrun}/  
	mv .//climb.genomics.cn/pub/10.5524/100001_101000/100051/data/data/${x}/BOLD/${task}_${myrun}/bold.nii.gz ./${x}/BOLD/${task}_${myrun}/bold.nii.gz 
      fi 
    done
  done
done 
