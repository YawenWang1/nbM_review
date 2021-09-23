#! /bin/bash
# This code is for registering UNI Image from FatNav-MP2RAGE sequence to MNI
# Colin27 space using antsAI and antsRegistrationSyN.sh
#  It took in 1h:6m:11s.
# Time the code
start_tme=$(date +%s)
# Start bold text
bold=$(tput bold)
# Turn off all attributes
normal=$(tput sgr0)
################################################################################
colin='/media/g/P01/Data/SPM_Anatomy/V22/colin27T1_seg.nii'
nbmPth='/media/g/P01/Data/SPM_Anatomy/V22/nbm_ROI/'
fatNavPth='/media/g/P04/Data/BIDS/'
TRFs='/media/g/P01/Data/SPM_Anatomy/V22/nbm_ROI/TRFs/'
Warps='/media/g/P01/Data/SPM_Anatomy/V22/nbm_ROI/Warps/'
mkdir ${Warps}
aryFolder=(1 \
"sub-07" \
"sub-08" \
"sub-10" \
"sub-11" \
"sub-12")
arySorder=(1 \
"sub-01" \
"sub-02" \
"sub-03" \
"sub-04" \
"sub-05")
strPst='/ses-002/anat/spm_bf_correction/'
################################################################################
# create non-csf and vessel brainmask
for subjd in {1..4}; do
  # Directory for the masks and brains
  currPth=${fatNavPth}${aryFolder[${subjd}]}${strPst}

  # Create brainmask without vessels
  fslmaths ${currPth}mINV2_brainmask.nii.gz -mul ${currPth}${aryFolder[${subjd}]}_non_brainmask_ivt.nii.gz ${nbmPth}${arySorder[${subjd}]}_brainmask_nov.nii.gz
  # Create brain masking out skull and vessels
  fslmaths ${currPth}mUNI_corrected.nii \
  -mul ${nbmPth}${arySorder[${subjd}]}_brainmask_nov.nii.gz \
  ${nbmPth}${arySorder[${subjd}]}_mUNI_corrected_brain.nii.gz

  # Estimate initial transformation matrix
  antsAI -d 3 \
  -m MI[ ${colin},${nbmPth}${arySorder[${subjd}]}_mUNI_corrected_brain.nii.gz,32,Regular,0.25] \
  -t AlignCentersOfMass \
  -o ${TRFs}${arySorder[${subjd}]}2Colin_init.mat
  # Register
  antsRegistrationSyN.sh -d 3 \
  -f ${colin} \
  -m ${nbmPth}${arySorder[${subjd}]}_mUNI_corrected_brain.nii.gz \
  -i ${TRFs}${arySorder[${subjd}]}2Colin_init.mat \
  -n 16 \
  -o ${Warps}${arySorder[${subjd}]}_to_colin_ \
  -t b \
  -e 13


done
echo " "
end_tme=$(date +%s)
nettme=$(expr ${end_tme} - ${start_tme})
echo "-----> It took in $((${nettme} / 3600))h:$((${nettme} % 3600 / 60))m:$((${nettme} % 60))s."
echo " "
#
# /opt/ANTs/bin/antsRegistration --verbose 1 --random-seed 13 --dimensionality 3 --float 0 --collapse-output-transforms 1 --output [ Sub-01_to_coline_tst_05_,Sub-01_to_coline_tst_05_Warped.nii.gz,Sub-01_to_coline_tst_05_InverseWarped.nii.gz ] --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [ 0.005,0.995 ] --initial-moving-transform Sub2Colin_init_01.mat --transform Rigid[ 0.1 ] --metric MI[ ../colin27T1_seg.nii,mUNI_corrected_brain_01.nii.gz,1,32,Regular,0.25 ] --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[ 0.1 ] --metric MI[ ../colin27T1_seg.nii,mUNI_corrected_brain_01.nii.gz,1,32,Regular,0.25 ] --convergence [ 1000x500x250x100,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform BSplineSyN[ 0.1,26,0,3 ] --metric CC[ ../colin27T1_seg.nii,mUNI_corrected_brain_01.nii.gz,1,4 ] --convergence [ 100x70x50x20,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
