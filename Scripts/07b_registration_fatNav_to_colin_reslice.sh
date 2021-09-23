#! /bin/bash
# This code is for registering UNI Image from FatNav sequence to MNI Colin space
# using antsAI and antsRegistrationSyN.sh
# Time the code
start_tme=$(date +%s)
# Start bold text
bold=$(tput bold)
# Turn off all attributes
normal=$(tput sgr0)
################################################################################
colin='/media/g/P01/Data/SPM_Anatomy/V22_anatoolbox/colin27T1_seg.nii'
nbmPth='/media/g/P01/Data/SPM_Anatomy/V22_anatoolbox/nbm_ROI/'
Warps='/media/g/P01/Data/SPM_Anatomy/V22_anatoolbox/nbm_ROI/Warps/'
SubROIs='/media/g/P01/Data/SPM_Anatomy/V22_anatoolbox/nbm_ROI/nbmPmaps/'
mkdir ${SubROIs}
RefnbM='/media/g/P01/Data/SPM_Anatomy/V22/Bforebrain_4.nii'
arySorder=(1 \
"sub-01" \
"sub-02" \
"sub-03" \
"sub-04" \
"sub-05")
################################################################################
# create non-csf and vessel brainmask
for subjd in {1..5}; do
  # Directory for the masks and brains
  currIn=${nbmPth}${arySorder[${subjd}]}_mUNI_corrected_brain.nii.gz
  trfMatrx=${Warps}${arySorder[${subjd}]}_to_colin_0GenericAffine.mat
  currOut=${SubROIs}${arySorder[${subjd}]}_nbm_pmap.nii.gz
  # Reslice
  $ANTSPATH/antsApplyTransforms \
  -d 3 \
  --float \
  -i ${RefnbM} \
  -r ${currIn} \
  -o ${currOut} \
  -n LanczosWindowedSinc \
  -t [${trfMatrx},1] \
  -v 1



done
echo " "
end_tme=$(date +%s)
nettme=$(expr ${end_tme} - ${start_tme})
echo "-----> It took in $((${nettme} / 3600))h:$((${nettme} % 3600 / 60))m:$((${nettme} % 60))s."
echo " "
