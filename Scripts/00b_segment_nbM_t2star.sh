#!/bin/bash

#-------------------------------------------------------------------------------
# This script is to segment nbM from ex vivo quantitative t2star image with the
# resolution of 0.2 mm isotropic.
currPth='/media/g/P01/Data/Alard'
# Generate folders that would be needed
for fldid in {0..5};
do
mkdir ${currPth}/op_0${fldid}
mkdir ${currPth}/L_smooth_0$(($fldid+1))
mkdir ${currPth}/R_smooth_0$(($fldid+1))
done

# Define the value of deltat and edgefraction value in an array. The first column
# is deltat and the second column is edgefraction
declare -A arr
arr[0,0]=0.2
arr[0,1]=0.7
arr[1,0]=0.2
arr[1,1]=0.7
arr[2,0]=0.3
arr[2,1]=0.7
arr[3,0]=0.1
arr[3,1]=0.7
arr[4,0]=0.1
arr[4,1]=0.8
arr[5,0]=0.1
arr[5,1]=0.6


# 3d 3danisosmooth
for fldid in {0..5};
do
cd ${currPth}/op_0${fldid}
3danisosmooth -deltat 0. -edgefraction 0.9 ../20160930_ExpT2DecGRE_masked_LIP_extract_00.nii.gz -prefix tst
3dAFNItoNIFTI -prefix ../s_deltat_pt$((arr[${fldid},0]*10))_edgefrac_pt$((arr[${fldid},1]*10)).nii SmoothAni+tlrc.BRIK.gz
cd ${currPth}
3dSeg -prefix R_smooth_0$(($fldid+1)) -mask Right_SI.nii -classes 'CSF ; GM ; WM' -bias_classes 'GM ; WM' -bias_fwhm 0.0 -anat s_deltat_pt$((arr[${fldid},0]*10))_edgefrac_pt$((arr[${fldid},1]*10)).nii
3dSeg -prefix L_smooth_0$(($fldid+1)) -mask Left_SI.nii -classes 'CSF ; GM ; WM' -bias_classes 'GM ; WM' -bias_fwhm 0.0 -anat s_deltat_pt$((arr[${fldid},0]*10))_edgefrac_pt$((arr[${fldid},1]*10)).nii
cd ${currPth}/L_smooth_0$(($fldid+1))
3dAFNItoNIFTI -prefix L_sm_t2star_op_0(($fldid+1))_Probabilisites.nii Posterior+tlrc.BRIK.gz
cd ${currPth}/R_smooth_0$(($fldid+1))
3dAFNItoNIFTI -prefix R_sm_t2star_op_0(($fldid+1))_Probabilisites.nii Posterior+tlrc.BRIK.gz
done

# segmentation without smoothing
3dSeg -prefix ${currPth}L_T2star_nbm_seg -mask ${currPth}Left_SI.nii -classes 'CSF ; GM ; WM' -bias_classes 'GM ; WM' -bias_fwhm 20 -anat ${currPth}20160930_ExpT2DecGRE_masked_LIP_extract_00.nii.gz
3dSeg -prefix ${currPth}R_T2star_nbm_seg -mask ${currPth}Right_SI.nii -classes 'CSF ; GM ; WM' -bias_classes 'GM ; WM' -bias_fwhm 20 -anat ${currPth}20160930_ExpT2DecGRE_masked_LIP_extract_00.nii.gz

# Convert afni files to nii.gz
3dAFNItoNIFTI -prefix ${currPth}L_T2star_nbm_seg/L_T2star_nbm_seg.nii.gz ${currPth}L_T2star_nbm_seg/Classes+tlrc.BRIK.gz
3dAFNItoNIFTI -prefix ${currPth}L_T2star_nbm_seg/L_T2star_nbm_Probabilities.nii.gz ${currPth}L_T2star_nbm_seg/Posterior+tlrc.BRIK.gz
3dAFNItoNIFTI -prefix ${currPth}R_T2star_nbm_seg/R_T2star_nbm_seg.nii.gz ${currPth}R_T2star_nbm_seg/Classes+tlrc.BRIK.gz
3dAFNItoNIFTI -prefix ${currPth}R_T2star_nbm_seg/R_T2star_nbm_Probabilities.nii.gz ${currPth}R_T2star_nbm_seg/Posterior+tlrc.BRIK.gz

# Examine the segmentation in ITK-SNAP
# When the deltat is 0.1 and edgefraction is 0.8, giving the best segmentation on the left
# side. The smoothing procedure didn't improve the segmentation of nbM on the
# right hemishpere
# Binarize the probabilisity maps
# Extract second cluster including nbM
fslroi ${currPth}L_smooth_05/L_sm_t2star_op_05_Probabilities.nii.gz ${currPth}L_smooth_05/L_sm_t2star_op_05_mask.nii.gz 1 1
fslroi ${currPth}R_T2star_nbm_seg/R_T2star_nbm_Probabilities.nii.gz ${currPth}R_T2star_nbm_seg/R_T2star_nbm_mask_01.nii.gz 1 1
# Binarize L hemishpere
fslmaths ${currPth}L_smooth_05/L_sm_t2star_op_05_mask.nii.gz -thr 0.5 -bin ${currPth}L_smooth_05/L_sm_t2star_op_05_mask_bin.nii.gz

# Binarize R hemishpere
fslmaths ${currPth}R_T2star_nbm_seg/R_T2star_nbm_mask.nii.gz -thr 0.01 -bin ${currPth}R_T2star_nbm_seg/R_T2star_nbm_mask_bin_01.nii.gz
fslmaths ${currPth}R_T2star_nbm_seg/R_T2star_nbm_mask.nii.gz -uthr 0.3 -bin ${currPth}R_T2star_nbm_seg/R_T2star_nbm_mask_bin_02.nii.gz
fslmaths ${currPth}R_T2star_nbm_seg/R_T2star_nbm_mask_bin_01.nii.gz -mul ${currPth}R_T2star_nbm_seg/R_T2star_nbm_mask_bin_02.nii.gz ${currPth}R_T2star_nbm_seg/R_T2star_nbm_mask_bin_03.nii.gz

# Adjust nbM mask in ITK-SNAP
