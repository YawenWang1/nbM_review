#! /bin/bash
################################################################################

# Mask brain mask from bias corrected INV2 Image
################################################################################
# Directory of sri's motion correction and distortion correction code
start_tme=$(date +%s)

################################################################################
# Simple formatting
# Start bold text
bold=$(tput bold)
# Turn off all attributes
normal=$(tput sgr0)

################################################################################


Pth='/media/h/P04/Data/BIDS/'
subject_id='sub-08'
session_id='ses-002'
TR=2.604
################################################################################
anatPth=${Pth}${subject_id}/${session_id}/anat/
anat_Pth=${Pth}${subject_id}/${session_id}/anat/spm_bf_correction
INV2=${anat_Pth}/mINV2_corrected.nii
INV1=${anat_Pth}/mINV1_corrected.nii
UNI=${anat_Pth}/mUNI_corrected.nii

# Decompress .nii.gz to .nii
# Get from sri
#-------------------------------------------------------------------------------
gunzip -c ${anat_Pth}/INV2_corrected.nii.gz > ${anat_Pth}/INV2_corrected.nii
gunzip -c ${anat_Pth}/INV1_corrected.nii.gz > ${anat_Pth}/INV1_corrected.nii
gunzip -c ${anat_Pth}/UNI_corrected.nii.gz > ${anat_Pth}/UNI_corrected.nii

# ################################################################################
# # Make mask
echo "-----> Creating Brain mask from INV2."
echo ""
/usr/share/fsl/5.0/bin/bet ${INV2} \
${INV2}_brain \
-f 0.20000000000000007 \
-g 0 \
-m \
-t

# Get non-brain mask
fslmaths \
${anat_Pth}/c3INV2_corrected.nii \
-add ${anat_Pth}/c4INV2_corrected.nii \
-add ${anat_Pth}/c5INV2_corrected.nii \
${anat_Pth}/sub-08_non_brainmask.nii.gz
# Manully adjust occipital lobe area to remove sinus
echo "Manully clean up around occipital lobe area"
# Get invertion of non_brain_Mask
fslmaths \
${anat_Pth}/sub-08_non_brainmask.nii.gz \
-sub 1 \
-mul -1 \
${anat_Pth}/sub-08_non_brainmask_ivt.nii.gz
# Calculate brains to fit to freesurfer
fslmaths \
${UNI} \
-mul ${anat_Pth}/mINV2_brainmask.nii.gz \
-mul ${anat_Pth}/sub-08_non_brainmask_ivt.nii.gz \
${anat_Pth}/sub-08_UNI_brain.nii.gz
# Get T1 brain
fslmaths \
${anat_Pth}/T1map_corrected.nii.gz \
-mul ${anat_Pth}/mINV2_brainmask.nii.gz \
${anat_Pth}/sub-08_T1_brain.nii.gz
# -mul ${anat_Pth}/sub-08_non_brainmask_ivt.nii.gz \
# ${anat_Pth}/sub-08_T1_brain.nii.gz


# Make a new folder for freesurfer files
fsPth=${Pth}${subject_id}/${session_id}/anat/fs
mkdir ${fsPth}
# Copy two files are needed for running freesurfer
cp ${anat_Pth}/sub-08_T1_brain.nii.gz ${fsPth}/sub-08_T1_brain.nii.gz
cp ${anat_Pth}/sub-08_UNI_brain.nii.gz ${fsPth}/sub-08_UNI_brain.nii.gz

cp /media/h/P04/Data/BIDS/sub-03/ses-002/anat/fs/expert.opts ${fsPth}/expert.opts



echo "-----> Finished to preprocess anatomical data."
echo " "
export SUBJECTS_DIR=${fsPth}
codePth='/media/h/P04/Data/BIDS/code/sub-08/'
subjid='sub-08'
${codePth}06c_run_MP2RAGE_Freesurfer_fixed.sh \
${subjid} \
${fsPth}/sub-08_UNI_brain.nii.gz \
${fsPth}/sub-08_T1_brain.nii.gz \
${fsPth}/expert.opts

echo " "
end_tme=$(date +%s)
nettme=$(expr $end_tme - $start_tme)
echo "-----> It took in $(($nettme / 3600))h:$(($nettme % 3600 / 60))m:$(($nettme % 60))s."
echo " "
