#!/bin/bash
#
dim=3 # image dimensionality
AP=""
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2  # controls multi-threading
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
f=$1 ; m=$2  ; bold=$3 ; studyid=$4 # fixed and moving image file names
if [[ ! -s $f ]] ; then echo no fixed $f ; exit; fi
if [[ ! -s $m ]] ; then echo no moving $m ;exit; fi
nm1=` basename $f | cut -d '.' -f 1 `
nm2=` basename $m | cut -d '.' -f 1 `
nm=${studyid}   # construct output prefix
reg=${AP}antsRegistration           # path to antsRegistration
echo affine $m $f outname is $nm
if [[ ! -s ${nm}_diff.nii.gz ]] ; then 
$reg -d $dim -r [ $f, $m ,1]  \
                        -m mattes[  $f, $m , 1 , 32 ] \
                         -t affine[ 0.1 ] \
                         -c [10000x1110x1100,1.e-7,10]  \
                        -s 2x1x0vox  -f 4x2x1  \
                        -m cc[  $f, $m , 1 , 4 ] \
                         -t syn[ 0.2,3,0 ] \
                         -c [40x40x0,0,10]  \
                        -s 2x1x0vox  -f 4x2x1 -l 1 -u 1 -z 1   \
                       -o [${nm},${nm}_diff.nii.gz]
  ${AP}antsApplyTransforms -d 3 -i $bold -r $f -t ${nm}1Warp.nii.gz -t ${nm}0GenericAffine.mat -o ${nm}.nii.gz -z 1 -e 3 
fi
