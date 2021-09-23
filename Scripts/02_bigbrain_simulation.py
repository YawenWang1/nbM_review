#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 10:06:32 2021
This script is for the simulation of rotation, translation and downsampling on 
BigBrain nbM ROI. The initial resolution is 0.1 mm isotropic

@author: yawen
"""
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
import nilearn
# from nilearn.image import resample_img
from scipy.ndimage import affine_transform
from scipy import ndimage
from scipy.ndimage.rotations import x_rotmat , y_rotmat, z_rotmat
import random
import pandas as pd
import seaborn as sns
import time
# =============================================================================
startime = time.time()
# =============================================================================
#------------------------------------------------------------------------------
# Define some functions


# def show_slices(slices):
#     """ Function to display row of image slices """
#     plt.figure(figsize=(2500/my_dpi,2500/my_dpi),dpi=my_dpi )
#     fig, axes = plt.subplots(1, len(slices))
#     for i, slice in enumerate(slices):
#         axes[i].imshow(slice.T, cmap="gray", origin="lower")

def is_even(x):
     if x % 2 == 0: 
          return True
     else:
          return False

#------------------------------------------------------------------------------

# =============================================================================
# Generate and save order of iteration as .npy file
random.seed(30)
order = []
for itid in range(Iteration_times):     
     curr_mot = random.randint(0,len(resolutions)-1)
     order.append(resolutions[curr_mot])

np.save(dataPth + 'Iterations_order.npy', order)      
# =============================================================================
my_dpi = 125
star, step, stop = 166, 10,1000  
reset_pt = np.arange(start=91,stop=1000,step=10)
# Get the original affine matrix 
L_nii_mask   = nib.load('/media/g/P01/Data/Bigbrain/Seg/L_BB_roi_nbm_seg/L_BB_roi_nbm_ROI_mask_bin_05.nii.gz')
R_nii_mask   = nib.load('/media/g/P01/Data/Bigbrain/Seg/R_BB_roi_nbm_seg/R_BB_roi_nbm_ROI_mask_bin_01.nii.gz')
# nii_brain  = nib.load('/media/g/P01/Data/Bigbrain/Seg/BB_roi_n.nii.gz')
dataPth = '/media/g/P01/Data/Bigbrain/Seg/tmp_mask_files/'
orders    = np.load('/media/g/P01/Data/Bigbrain/Seg/tmp_mask_files/Iterations_order.npy')
ary_L_nii_mask = L_nii_mask.dataobj[:]
ary_R_nii_mask = R_nii_mask.dataobj[:]

old_resolution=0.1
print("Original Resolution")
print(old_resolution)
resolutions = np.arange(start=0.1,stop=3.1,step=0.1)
resolutions = np.ndarray.round(resolutions,decimals=2)

nme_mask_prefix = 'mask_iteration_'
nme_brain_prefix = 'brain_iteration_'

del L_nii_mask
del R_nii_mask
# get_ipython().magic('reset -sf')           
# first stopping point :166   

lst_Mot_tmp = []
Rot_mask_tmp = []

for itid in range(star,stop):
     
     if itid in reset_pt:
          Rot_mask_tmp_vox, Rot_mask_tmp_vo1 = [],[]      

          for i in range(len(Rot_mask_tmp)):
               Rot_mask_tmp_vox.append(np.count_nonzero(Rot_mask_tmp[i]))
               Rot_mask_tmp_vo1.append(np.count_nonzero(Rot_mask_tmp[i])*(0.1**3))
      
          np.save(dataPth + 'Iter_' +str(itid-10) + '_' + str(itid-1) + '_num_vox.npy', Rot_mask_tmp_vox)          
          np.save(dataPth + 'Iter_' +str(itid-10) + '_' + str(itid-1) + '_num_vol.npy', Rot_mask_tmp_vo1)   
          
          # Delete big variable
          del lst_Mot_tmp
          del Rot_mask_tmp_vox
          del Rot_mask_tmp_vo1
          
          Num_voxels = np.zeros([len(Rot_mask_tmp),len(resolutions)])
          Vol_size = np.zeros([len(Rot_mask_tmp),len(resolutions)])
          lst_num_vox = []
          lst_vol_sze = []
          lst_lr      = []
          x_axis      = []
          lst_ds      = []
          lst_old_vol = []
          lst_old_num = []
          for i in range(len(Rot_mask_tmp)):
               count = -1
          
               if is_even(np.divmod(i,6)[0]):
          
                 lr     = 'Left'

               else:

                 lr     = 'Right'
          
               for new_resolution in resolutions:
                   # lst_mot.append(lr_mot)
                   lst_lr.append(lr)
                   count += 1
                   print("New Resolution")
                   print(new_resolution)
                   downsample_factor = new_resolution / old_resolution
                   print("Downsampling Factor")
                   print(downsample_factor)
               
                   old_size=np.asarray((Rot_mask_tmp[i].shape[0],Rot_mask_tmp[i].shape[1],Rot_mask_tmp[i].shape[2]),dtype=int)
                   print("-------------------------------------------------")
                   print("Old Dimension")
                   print(tuple(old_size))
                   new_size=tuple(old_size // downsample_factor)
                   print("New Dimension")
                   print(tuple(new_size))
               
                   print("Downsampling Data")
                   # new_data_img=ndimage.zoom(input=Rot_mask[i],zoom=1/downsample_factor,order=5)
                   # np.save(str(data_fname+"_ds_"+str(round(downsample_factor,3))+"x.npy"),new_data_img)
                   new_mask_img=ndimage.zoom(input=Rot_mask_tmp[i],zoom=1/downsample_factor,order=0,mode='nearest')
                   lst_ds.append(new_mask_img)
                   # np.save(str(mask_fname+"_ds_"+str(round(downsample_factor,3))+"x.npy"),new_mask_img)
               
                   print("-------------------------------------------------")
                   print("ROI voxels")
                   old_mask_voxels = np.count_nonzero(Rot_mask_tmp[i])
                   new_mask_voxels = np.count_nonzero(new_mask_img)
                   print(old_mask_voxels,new_mask_voxels)
                   Num_voxels[i,count] = new_mask_voxels
                   lst_num_vox.append(new_mask_voxels)
                   print("Voxel volume")
                   old_voxel_volume = 0.1*0.1*0.1 
                   new_voxel_volume = old_voxel_volume * (downsample_factor **3)
                   print(old_voxel_volume,new_voxel_volume)
               
                   print("ROI Volume")
                   print(old_mask_voxels*old_voxel_volume,new_mask_voxels*new_voxel_volume)
                   Vol_size[i,count] = new_mask_voxels*new_voxel_volume
                   lst_vol_sze.append(new_mask_voxels*new_voxel_volume)
                   lst_old_vol.append(old_mask_voxels*old_voxel_volume)
                   lst_old_num.append(old_mask_voxels)
                   x_axis.append(new_resolution)
                   print("-------------------------------------------------")
                   
          df_tmp = pd.DataFrame(data={'Hemisphere': lst_lr,
                        'Volume Size': lst_vol_sze,
                        'Num of Voxels': lst_num_vox,
                        'Resolution': x_axis,
                        'Orig_num_voxel': lst_old_num,
                        'Orig_vol_size': lst_old_vol})         
          df_tmp.to_excel(dataPth + 'Iter_' +str(itid-10) + '_' + str(itid-1) + '_downsampling_df.xlsx')
          
          del Rot_mask_tmp
     
       
          plt.figure(figsize=(2500/my_dpi,2500/my_dpi),dpi=my_dpi )
          ax4 = sns.barplot(data=df_tmp,x="Resolution",y="Volume Size",hue="Hemisphere",palette="mako")
          # ax6.axhline(np.count_nonzero(L_nii_mask_ary)*(0.1**3),linestyle='-.', color='black',label='LH_BigBrain_NBM')
          # ax6.axhline(np.count_nonzero(R_nii_mask_ary)*(0.1**3),linestyle='-.',color='green',label='RH_BigBrain_NBM')
          ax4.axes.set_title("Volume size: downsampling of nbM with Rotation and Translation",fontsize=30)
          ax4.set_xlabel("Spatial resolution",fontsize=20)
          ax4.set_ylabel("Volume Size",fontsize=20)
          ax4.tick_params(labelsize=13)
          plt.savefig(dataPth + 'Iter_' +str(itid-10) + '_' + str(itid-1)  + '_downsampling_vol.png',dpi=300)
          plt.savefig(dataPth + 'Iter_' +str(itid-10) + '_' + str(itid-1)  + '_downsampling_vol.pdf',dpi=300)
          
          plt.close('all')
          
          plt.figure(figsize=(2500/my_dpi,2500/my_dpi),dpi=my_dpi )
          ax5 = sns.barplot(data=df_tmp,x="Resolution",y="Num of Voxels",hue="Hemisphere",palette="mako")
          # ax6.axhline(np.count_nonzero(L_nii_mask_ary)*(0.1**3),linestyle='-.', color='black',label='LH_BigBrain_NBM')
          # ax6.axhline(np.count_nonzero(R_nii_mask_ary)*(0.1**3),linestyle='-.',color='green',label='RH_BigBrain_NBM')
          ax5.axes.set_title("Number of Voxel: downsampling of nbM with Rotation and Translation",fontsize=30)
          ax5.set_xlabel("Spatial resolution",fontsize=20)
          ax5.set_ylabel("Volume Size",fontsize=20)
          ax5.tick_params(labelsize=13)
          plt.savefig(dataPth + 'Iter_' +str(itid-10) + '_' + str(itid-1) + '_downsampling_vol_01_30.png',dpi=300)
          plt.savefig(dataPth + 'Iter_' +str(itid-10) + '_' + str(itid-1) + '_downsampling_vol_01_30.pdf',dpi=300)
          
          plt.close('all')
          
          
          df_tmp_partial = df_tmp.loc[df_tmp['Resolution']>0.4]
          plt.figure(figsize=(2500/my_dpi,2500/my_dpi),dpi=my_dpi )
          ax6 = sns.barplot(data=df_tmp_partial,x="Resolution",y="Num of Voxels",hue="Hemisphere",palette="mako")
          # ax6.axhline(np.count_nonzero(L_nii_mask_ary)*(0.1**3),linestyle='-.', color='black',label='LH_BigBrain_NBM')
          # ax6.axhline(np.count_nonzero(R_nii_mask_ary)*(0.1**3),linestyle='-.',color='green',label='RH_BigBrain_NBM')
          ax6.axes.set_title("Number of Voxel: downsampling of nbM with Rotation and Translation",fontsize=30)
          ax6.set_xlabel("Spatial resolution",fontsize=20)
          ax6.set_ylabel("Volume Size",fontsize=20)
          ax6.tick_params(labelsize=13)
          plt.savefig(dataPth + 'Iter_' +str(itid-10) + '_' + str(itid-1) + '_downsampling_vol_05_30.png',dpi=300)
          plt.savefig(dataPth + 'Iter_' +str(itid-10) + '_' + str(itid-1) + '_downsampling_vol_05_30.pdf',dpi=300)
          
          plt.close('all')
          
          del df_tmp
          del df_tmp_partial
          del ax4
          del ax5 
          del ax6
          del Num_voxels 
          del Vol_size 
          del lst_num_vox 
          del lst_vol_sze 
          del lst_lr     
          del x_axis     
          del lst_ds      
          del lst_old_vol
          del lst_old_num 
          
          lst_Mot_tmp = []
          Rot_mask_tmp = []
                         


     curr_mot = orders[itid]
     Mot = np.array([[curr_mot,curr_mot,curr_mot],
                  [curr_mot,-curr_mot,curr_mot],
                  [curr_mot,-curr_mot,-curr_mot],
                  [-curr_mot,-curr_mot,-curr_mot],
                  [-curr_mot,curr_mot,-curr_mot],
                  [-curr_mot,curr_mot,curr_mot]])
     Rot_mat = np.ndarray.round(np.deg2rad(Mot),decimals=6)
  
     for lrid in range(2):
        if lrid == 0:
            curr_mask =  ary_L_nii_mask
        else:
            curr_mask =  ary_R_nii_mask
            
        for mid in range(len(Mot)):
    
            rot_matrix =  np.matmul(np.matmul(x_rotmat(Rot_mat[mid,0]),y_rotmat(Rot_mat[mid,1])),z_rotmat(Rot_mat[mid,2]))
            translation = np.round((Mot[mid,:]/0.1))
            curr_rot = affine_transform(curr_mask, rot_matrix,translation,mode='nearest')
            Rot_mask_tmp.append(curr_rot)
            np.save(dataPth + nme_mask_prefix + str(itid) +
                    '_rot_' + '_'.join(x for x in Rot_mat[mid,:].astype(str)) 
                    + '_trans_' +'_'.join(x for x in translation.astype(str)) 
                    +'.npy',curr_rot)
            lst_Mot_tmp.append(np.hstack([Rot_mat[mid,:],translation]))    
            
     

     
endtime = time.time()

durtime = endtime - startime 



print(f'This code takes: {durtime} seconds to run')





