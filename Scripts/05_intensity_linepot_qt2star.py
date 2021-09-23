#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 07:55:28 2021

@author: yawen
"""
import os
import glob
import nibabel as nib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from nibabel import load, Nifti1Image, save
from matplotlib.ticker import LinearLocator
from matplotlib import cm
from sklearn import preprocessing
# Define some parameters
FigPth = '/media/g/P01/Review/Figures/' 
extnii = load('/media/g/P01/Data/Alard/20160930_ExpT2DecGRE_masked_LIP_extract_00.nii.gz')
nbmnii = load('/media/g/P01/Data/Alard/L_smooth_05/L_sm_t2star_op_05_mask_nbm.nii.gz')
ary_extnii, ary_nbmnii = extnii.get_data(), nbmnii.get_data()

# Extract the slice that would be cropped in both Brain image and nbM ROI image

Image_s46, nbM_s46 = ary_extnii[46,:,:], ary_nbmnii[46,:,:]

idx_nonzero = np.nonzero(nbM_s46)

row_idx,row_idx_c = np.unique(idx_nonzero[0], return_counts=True)

clm_idx,clm_idx_c = np.unique(idx_nonzero[1], return_counts=True)


#%% Plot the images on the black background
norm_img_s46 = preprocessing.normalize(Image_s46)

# Crop a small piece of image from the whole brain
# Set up the range. 
ext = 10
y_range = [clm_idx[0] - ext, clm_idx[-1] + ext]
x_range = [row_idx[0] - ext, row_idx[-1] + ext]

norm_img_s46_crop = norm_img_s46[x_range[0]:x_range[1],y_range[0]:y_range[1]]


fig = plt.figure(figsize=(ary_extnii.shape[1]/my_dpi,ary_extnii.shape[2]/my_dpi),dpi=my_dpi )
ax = fig.add_axes([0, 0, 1, 1])
# Hide spines, ticks, etc.
ax.axis('off')
boundaries_01 = np.linspace(10,57,10)
# ax0.imshow(norm_img_s46_crop,cmap='gray',vmin=0.00660,vmax=0.05)
# set up fontsize
ftsz = 16
my_dpi = 125

rc_hist={'font.size': ftsz,
 'axes.labelsize': ftsz,
 'axes.titlesize': ftsz,
 'xtick.labelsize': ftsz,
 'ytick.labelsize': ftsz,
 'xtick.color':'white', 'ytick.color':'white',
 'legend.fontsize': 7.0,
 'axes.linewidth': 1.25,
 'grid.linewidth': 1.0,
 'lines.linewidth': 1.5,
 'lines.markersize': 6.0,
 'patch.linewidth': 1.0,
 'xtick.major.width': 1,
 'ytick.major.width': 1,
 'xtick.minor.width': 1.0,
 'ytick.minor.width': 1.0,
 'xtick.major.size': 1.0,
 'ytick.major.size': 1.0,
 'xtick.minor.size': 1.0,
 'ytick.minor.size': 1.0,
 'legend.title_fontsize': 8.0,
 'font.family':'Arial'}

fig, [ax0,ax1] = plt.subplots(1,2,figsize = [8.06,6], gridspec_kw={'width_ratios': [1.3, 1]}, sharex = False, sharey= True, dpi=my_dpi)
plt.rcParams.update(**rc_hist)
fig.patch.set_facecolor('xkcd:black')
plt.subplots_adjust(left =0.1, right = 0.93,top = 0.975, wspace = 0.134, bottom=0.03 )
ax0_cp = ax0.twinx()  
ax0.plot(np.mean(norm_img_s46_crop[10:-10,:],axis=0),color='slateblue',linewidth=4.0)
ax0.annotate('Horizontal: L to R', (30, 0.045),fontsize=ftsz)
ax0.set_ylim(0,0.05)
ax0.set_ylabel('Normalized intensity', color='slateblue', fontsize=ftsz)
counts, bins = np.histogram(idx_nonzero[1]-clm_idx[0]+10, bins=np.arange(norm_img_s46_crop.shape[1]))
factor = 0.01
#1 ax0_cp.hist(bins[:-1], bins, weights=factor*counts,color='grey', alpha=0.3)
ax0_cp.hist(bins[:-1], bins, weights=counts, color='grey', alpha=0.3)
ax0_cp.set_ylim(0,28.5)
# ax0.fill_between(np.arange(10,56),0.05, color='grey', alpha=0.3)
ax1_cp = ax1.twinx()  
ax1.plot(np.mean(norm_img_s46_crop[:,10:-10],axis=1),color='slateblue',linewidth=4.0)
ax1.set_ylim(0,0.05)
ax1.annotate('Vertical: I to S', (10, 0.045),fontsize=ftsz)
counts_01, bins_01 = np.histogram(idx_nonzero[0]-row_idx[0]+10, bins=np.arange(norm_img_s46_crop.shape[0]))
ax1_cp.hist(bins_01[:-1], bins_01, weights=counts_01, color='grey', alpha=0.3)
ax1_cp.set_ylabel('Histogram of the number of nbM voxels', color='grey', fontsize=ftsz)
ax1_cp.set_ylim(0,28.5)
plt.savefig(FigPth + 'nbM review figure 5C k bg 01.pdf',dpi=300)

plt.savefig(FigPth + 'nbM review figure 5C k bg 01.svg',dpi=300)

plt.savefig(FigPth + 'nbM review figure 5C k bg 01.png',dpi=300)

# =============================================================================
# # Plot the images on the white background
# rc_hist={'font.size': ftsz,
#  'axes.labelsize': ftsz,
#  'axes.titlesize': ftsz,
#  'xtick.labelsize': ftsz,
#  'ytick.labelsize': ftsz,
#  'legend.fontsize': 7.0,
#  'axes.linewidth': 1.25,
#  'grid.linewidth': 1.0,
#  'lines.linewidth': 1.5,
#  'lines.markersize': 6.0,
#  'patch.linewidth': 1.0,
#  'xtick.major.width': 1,
#  'ytick.major.width': 1,
#  'xtick.minor.width': 1.0,
#  'ytick.minor.width': 1.0,
#  'xtick.major.size': 1.0,
#  'ytick.major.size': 1.0,
#  'xtick.minor.size': 1.0,
#  'ytick.minor.size': 1.0,
#  'legend.title_fontsize': 8.0,
#  'font.family':'Arial'}


# fig, [ax0,ax1] = plt.subplots(1,2,figsize = [8.06,6], gridspec_kw={'width_ratios': [1.3, 1]}, sharex = False, sharey= True, dpi=my_dpi)
# plt.rcParams.update(**rc_hist)

# plt.subplots_adjust(left =0.1, right = 0.93,top = 0.975, wspace = 0.134, bottom=0.03 )
# ax0_cp = ax0.twinx()  
# ax0.plot(np.mean(norm_img_s46_crop[10:-10,:],axis=0),color='black',linewidth=4.0)
# ax0.annotate('Horizontal: L to R', (30, 0.045),fontsize=ftsz)
# ax0.set_ylim(0,0.05)
# ax0.set_ylabel('Normalized intensity', color='black', fontsize=ftsz)
# counts, bins = np.histogram(idx_nonzero[1]-clm_idx[0]+10, bins=np.arange(norm_img_s46_crop.shape[1]))
# factor = 0.01
# #1 ax0_cp.hist(bins[:-1], bins, weights=factor*counts,color='grey', alpha=0.3)
# ax0_cp.hist(bins[:-1], bins, weights=counts, color='grey', alpha=0.3)
# # ax0.fill_between(np.arange(10,56),0.05, color='grey', alpha=0.3)
# ax1_cp = ax1.twinx()  
# ax1.plot(np.mean(norm_img_s46_crop[:,10:-10],axis=1),color='black',linewidth=4.0)
# ax1.set_ylim(0,0.05)
# ax1.annotate('Vertical: I to S', (10, 0.045),fontsize=ftsz)
# counts_01, bins_01 = np.histogram(idx_nonzero[0]-row_idx[0]+10, bins=np.arange(norm_img_s46_crop.shape[0]))
# ax1_cp.hist(bins_01[:-1], bins_01, weights=counts_01, color='grey', alpha=0.3)
# ax1_cp.set_ylabel('Histogram of the number of nbM voxels', color='grey', fontsize=ftsz)



# plt.savefig(FigPth + 'nbM review figure 5C.pdf',dpi=300)

# plt.savefig(FigPth + 'nbM review figure 5C.svg',dpi=300)

# plt.savefig(FigPth + 'nbM review figure 5C.png',dpi=300)

# =============================================================================
