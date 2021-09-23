#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 13:17:42 2021
This script is for extracting left and right nbm from Julich brain version 2.2

@author: yawen
"""

import nibabel as nib
import numpy as np
# ----------------------------------------------------------------------------------
# Julich brain MNI Colin27 space

mpmPth='/media/g/P01/Data/SPM_Anatomy/V22_Jubrain/MPM/'
RH =  (mpmPth + 
     'JulichBrain_MPMAtlas_r_N10_nlin2Stdcolin27_publicDOI_862e709c1258b99106d936bbaccf0bdb.nii.gz')
LH =  (mpmPth + 
     'JulichBrain_MPMAtlas_l_N10_nlin2Stdcolin27_publicDOI_4514e5af7a15c96e47a049a98444ed7a.nii.gz')

basename = 'V22_jubrain_colin27'

rh_nii = nib.load(RH)

lh_nii = nib.load(LH)

ary_rh = rh_nii.get_data()

ary_lh = lh_nii.get_data()

nbm_value = 12

ary_rh[np.where(ary_rh!=nbm_value)], ary_rh[np.where(ary_rh==nbm_value)] =0 , 1

ary_lh[np.where(ary_lh!=nbm_value)], ary_lh[np.where(ary_lh==nbm_value)] =0 , 1


img_rnbm = nib.Nifti1Image(ary_rh, header=rh_nii.header, affine=rh_nii.affine)
nib.save(img_rnbm, mpmPth + basename + "_r_ch4.nii")



img_lnbm = nib.Nifti1Image(ary_lh, header=lh_nii.header, affine=lh_nii.affine)
nib.save(img_lnbm, mpmPth + basename + "_l_ch4.nii")


ary_ch4 = np.add(ary_lh,ary_rh)

img_nbm = nib.Nifti1Image(ary_ch4, header=lh_nii.header, affine=lh_nii.affine)
nib.save(img_nbm, mpmPth  + "V22_bl_jubrain_colin27_ch4.nii")

#------------------------------------------------------------------------------
# JUlich brain MNI ICBM152 space
RH =  (mpmPth + 
     'JulichBrain_MPMAtlas_r_N10_nlin2Stdicbm152asym2009c_publicDOI_14622b49a715338ce96e96611d395646.nii.gz')
LH =  (mpmPth + 
     'JulichBrain_MPMAtlas_l_N10_nlin2Stdicbm152asym2009c_publicDOI_3f6407380a69007a54f5e13f3c1ba2e6.nii.gz')

basename_icbm = 'V22_jubrain_icbm152'


rh_nii = nib.load(RH)

lh_nii = nib.load(LH)

ary_rh = rh_nii.get_data()

ary_lh = lh_nii.get_data()


ary_rh[np.where(ary_rh!=nbm_value)], ary_rh[np.where(ary_rh==nbm_value)] =0 , 1

ary_lh[np.where(ary_lh!=nbm_value)], ary_lh[np.where(ary_lh==nbm_value)] =0 , 1


img_rnbm = nib.Nifti1Image(ary_rh, header=rh_nii.header, affine=rh_nii.affine)
nib.save(img_rnbm, mpmPth + basename_icbm + "_r_ch4.nii")



img_lnbm = nib.Nifti1Image(ary_lh, header=lh_nii.header, affine=lh_nii.affine)
nib.save(img_lnbm, mpmPth + basename_icbm + "_l_ch4.nii")

ary_ch4 = np.add(ary_lh,ary_rh)

img_nbm = nib.Nifti1Image(ary_ch4, header=lh_nii.header, affine=lh_nii.affine)
nib.save(img_nbm, mpmPth + "V22_bl_jubrain_icbm152_ch4.nii")
