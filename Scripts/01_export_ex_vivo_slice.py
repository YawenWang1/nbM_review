#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 10:06:32 2021
This script is for extract one slice from qT2star, MGH brain, ding altas and  
to make figure 5, and supplementary figure 3 and 4.
@author: yawen
"""
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib

# =============================================================================
# % MGH brain
brain = '/media/g/P01/Review/Videos/synthesized_FLASH25_in_MNI.nii.gz'

nii = nib.load(brain)


my_dpi = 125

plt.figure(figsize=(nii.dataobj.shape[0]/my_dpi,nii.dataobj.shape[2]/my_dpi),dpi=my_dpi )
plt.imshow(nii.dataobj[:,889,:],cmap='gray',vmin=8.75,vmax=24)
plt.savefig('/media/g/P01/Review/Figures/synthesized_FLASH25_in_MNI_coronal_slice_'+str(889)+'.pdf',dpi=2000)
plt.savefig('/media/g/P01/Review/Figures/synthesized_FLASH25_in_MNI_coronal_slice_'+str(889)+'.png',dpi=2000)
# =============================================================================


# =============================================================================
# quantitative T2star

T2star = '/media/g/P01/Data/Alard/20160930_ExpT2DecGRE_masked_LIP_extract_00.nii.gz'
T2starniii = nib.load(T2star)
sl = 46

# plot
plt.figure(figsize=(T2starniii.dataobj.shape[0]/my_dpi,T2starniii.dataobj.shape[2]/my_dpi),dpi=my_dpi )
plt.imshow(T2starniii.dataobj[sl,:,:],cmap='gray',vmin=0.007,vmax=0.0356)
plt.savefig('/media/g/P01/Review/Figures/20160930_ExpT2DecGRE_masked_LIP_'+str(sl)+'.pdf',dpi=2000)
plt.savefig('/media/g/P01/Review/Figures/20160930_ExpT2DecGRE_masked_LIP_'+str(sl)+'.png',dpi=2000)
plt.close('all')
      
# =============================================================================
     
# =============================================================================
# Ding atlas

dingPth = '/media/h/Atlas/Ding/'

ding_brain = ['flash20_rot_crop',
          'flash40_rot_crop',
          'flash60_rot_crop',
          'flash80_rot_crop',
          'PD_rot_crop',
          'T2star_rot_crop',
          'T1_rot_crop']

# define the coronal slice that being exported and resolution of the screen
sl, my_dpi = 163, 125
aryVmax = [300.0, 270.0, 260.0, 265.0, 1000, 260, 220]
aryVmin = [74.5, 74.5, 74.5, 74.5, 74.5, 30, 30]
counter = -1
for i in ding_brain:
     counter += 1
     nii = nib.load(dingPth + i + '.nii.gz')
     plt.figure(figsize=(nii.dataobj.shape[0]/my_dpi,nii.dataobj.shape[2]/my_dpi),dpi=my_dpi )
     plt.imshow(nii.dataobj[:,sl,:],cmap='gray',vmin=aryVmin[counter],vmax=aryVmax[counter])
     # plt.savefig('/media/g/P01/Review/Figures/ding_'+ i + '_' + str(sl)+'.pdf',dpi=2000)
     # plt.savefig('/media/g/P01/Review/Figures/ding_'+ i + '_' + str(sl)+'.png',dpi=2000)
     plt.savefig('/media/g/P01/Review/Figures/ding_'+ i + '_' + str(sl)+'_dpi_127.pdf',dpi=127)
     plt.savefig('/media/g/P01/Review/Figures/ding_'+ i + '_' + str(sl)+'_dpi_127.png',dpi=127)

plt.close('all')

# =============================================================================
      

     