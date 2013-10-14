#!/bin/bash
#
dim=3 # image dimensionality
AP=""
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2  # controls multi-threading
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
f=$1 ; m=$2  ; bold1=$3  ;  bold2=$4 ; studyid=$5 # fixed and moving image file names
if [[ ! -s $f ]] ; then echo no fixed $f ; exit; fi
if [[ ! -s $m ]] ; then echo no moving $m ;exit; fi
nm1=` basename $f | cut -d '.' -f 1 `
nm2=` basename $m | cut -d '.' -f 1 `
nm=bold2to1_${studyid}   # construct output prefix
reg=${AP}antsRegistration           # path to antsRegistration
echo affine $m $f outname is $nm
$reg -d $dim -r [ $f, $m ,1]  \
                        -m mattes[  $f, $m , 1 , 32, regular, 0.25 ] \
                         -t rigid[ 0.2 ] \
                         -c 10000x1110x100  \
                        -s 2x1x0vox  -f 4x2x1 -l 1 -u 1 -z 1  \
                       -o [${nm},${nm}_aff.nii.gz]
#
${AP}antsApplyTransforms -d $dim -i $bold2 -r $f -n NearestNeighbor -t ${nm}0GenericAffine.mat -o ${nm}_warped.nii.gz -z 1
MeasureImageSimilarity $dim 0 $bold1 ${nm}_warped.nii.gz ${nm}_log.txt
MeasureImageSimilarity $dim 1 $bold1 ${nm}_warped.nii.gz ${nm}_log.txt
ThresholdImage 3 ${nm}_warped.nii.gz ${nm}_evalb2.nii.gz 0.000001 9
ThresholdImage 3 $bold1 ${nm}_evalb.nii.gz 0.000001 9
ImageMath 3 ${nm}_eval DiceAndMinDistSum ${nm}_evalb2.nii.gz ${nm}_evalb.nii.gz ${nm}_dist.nii.gz\