#!/bin/bash

#-------------------------------------------------------------------------------
#This script is to segment nbM from BigBrain dataset using afni tools

start_tme=$(date +%s)




bbrainPth='/media/g/P01/Data/Bigbrain/'
# Crop a small region to reduce computational time. The small region includes
# the whole nbM and the neighboring structures.

# fslroi ${bbrainPth}full8_100um_2009b_sym_extract.nii.gz ${bbrainPth}Seg/BB_roi_n.nii.gz 500 900 900 500 10 140

# Extract the brain images on both hemishphere based on the masks on the both Hemisphere
fslmaths ${bbrainPth}Seg/BB_roi_n.nii.gz -mul ${bbrainPth}Seg/L_BB_roi_04.nii.gz ${bbrainPth}Seg/L_BB_snbm.nii.gz
fslmaths ${bbrainPth}Seg/BB_roi_n.nii.gz -mul ${bbrainPth}Seg/R_BB_roi_01.nii.gz ${bbrainPth}Seg/R_BB_snbm.nii.gz

# Segment 'CSF ; GM ; WM'
3dSeg -prefix ${bbrainPth}Seg/L_BB_roi_nbm_seg -mask ${bbrainPth}Seg/L_BB_roi_03.nii.gz -classes 'CSF ; GM ; WM' -bias_classes 'GM ; WM' -bias_fwhm 20 -anat ${bbrainPth}Seg/L_BB_snbm.nii.gz
3dSeg -prefix ${bbrainPth}Seg/R_BB_roi_nbm_seg -mask ${bbrainPth}Seg/R_BB_roi_01.nii.gz -classes 'CSF ; GM ; WM' -bias_classes 'GM ; WM' -bias_fwhm 20 -anat ${bbrainPth}Seg/R_BB_snbm.nii.gz

# Convert afni files to nii.gz
3dAFNItoNIFTI -prefix ${bbrainPth}Seg/L_BB_roi_nbm_seg/L_BB_roi_nbm_Seg.nii.gz ${bbrainPth}Seg/L_BB_roi_nbm_seg/Classes+tlrc.BRIK.gz
3dAFNItoNIFTI -prefix ${bbrainPth}Seg/L_BB_roi_nbm_seg/L_BB_roi_nbm_Probabilities.nii.gz ${bbrainPth}Seg/L_BB_roi_nbm_seg/Posterior+tlrc.BRIK.gz
3dAFNItoNIFTI -prefix ${bbrainPth}Seg/R_BB_roi_nbm_seg/R_BB_roi_nbm_Seg.nii.gz ${bbrainPth}Seg/R_BB_roi_nbm_seg/Classes+tlrc.BRIK.gz
3dAFNItoNIFTI -prefix ${bbrainPth}Seg/R_BB_roi_nbm_seg/R_BB_roi_nbm_Probabilities.nii.gz ${bbrainPth}Seg/R_BB_roi_nbm_seg/Posterior+tlrc.BRIK.gz

# Extract the cluster including nbM
fslroi ${bbrainPth}Seg/L_BB_roi_nbm_seg/L_BB_roi_nbm_Probabilities.nii.gz ${bbrainPth}Seg/L_BB_roi_nbm_seg/L_BB_roi_nbm_ROI_mask.nii.gz 0 1
fslroi ${bbrainPth}Seg/R_BB_roi_nbm_seg/R_BB_roi_nbm_Probabilities.nii.gz ${bbrainPth}Seg/R_BB_roi_nbm_seg/R_BB_roi_nbm_ROI_mask.nii.gz 0 1
# Binarize them
fslmaths ${bbrainPth}Seg/L_BB_roi_nbm_seg/L_BB_roi_nbm_ROI_mask.nii.gz -thr 0.6 -bin ${bbrainPth}Seg/L_BB_roi_nbm_seg/L_BB_roi_nbm_ROI_mask_bin.nii.gz
fslmaths ${bbrainPth}Seg/R_BB_roi_nbm_seg/R_BB_roi_nbm_ROI_mask.nii.gz -thr 0.6 -bin ${bbrainPth}Seg/R_BB_roi_nbm_seg/R_BB_roi_nbm_ROI_mask_bin.nii.gz
# fslmaths BB_roi_nbm_ROI_mask.nii.gz -thr 0.5 -bin BB_roi_nbm_ROI_mask.nii.gz
# itksnap -g BB_roi_nbm.nii.gz -s BB_roi_nbm_Seg.nii.gz
#
#
# BB_roi_nbm_thr.nii.gz -> fslmaths  BB_roi_nbm.nii.gz -thr 3e4  BB_roi_nbm_thr.nii.gz
# BB_roi_nbm_mask.nii.gz -> fslmaths  BB_roi_nbm.nii.gz -bin  BB_roi_nbm_mask.nii.gz


end_tme=$(date +%s)
nettme=$(expr $end_tme - $start_tme)
echo " "
echo "-----> the automatic segmentation of nbM took in $(($nettme / 3600))h:$(($nettme % 3600 / 60))m:$(($nettme % 60))s."
echo " "
