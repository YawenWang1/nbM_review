#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 14:48:44 2021
This scirpts is for thresholding probabilistic maps from spm anatomy toolbox version 2.2
on individual level
@author: yawen
"""
import os
import glob
import nibabel as nib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
# =============================================================================
# *** Define functions
def fncLoadNii(strPathIn):
    """Load nii files."""
    print(('---------Loading: ' + strPathIn))
    # Load nii file (this doesn't load the data into memory yet):
    niiTmp = nib.load(strPathIn)
    # Load data into array:
    aryTmp = niiTmp.get_fdata()
    # Get headers:
    hdrTmp = niiTmp.header
    # Get 'affine':
    niiAff = niiTmp.affine
    # Output nii data as numpy array and header:
    return aryTmp, hdrTmp, niiAff

pMapsPth = '/media/g/P01/Data/SPM_Anatomy/V22_anatoolbox/nbm_ROI/nbmPmaps/'
figPth = '/media/g/P01/Review/Figures/'
os.chdir(pMapsPth)

AllFiles = glob.glob('*.nii.gz')
upper, lower = 1, 0


Voxel_count = []
Vox_c1d = []
subjs_nme = ['sub-01','sub-02','sub-03','sub-04','sub-05']
# threshold = [0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7]

threshold = [0.3,0.4,0.5,0.6,0.7]
Thre = []
Subjs = []
Volume=[]
DataType=[]
for fileid in range(0,len(AllFiles)):
     
     nii, hdrnii, affnii = fncLoadNii(AllFiles[fileid])
     Num = []
     for th in threshold:
          num = np.count_nonzero(np.where(nii>th,upper,lower))
          Num.extend([num])
          Vox_c1d.extend([num])
          Volume.extend([num*(0.8**3)])
          Thre.extend([th])
          Subjs.extend([AllFiles[fileid].split('_')[0]])
          DataType.extend(['in vivo 7T'])
          
     Voxel_count.append(Num)
     
pmaps_dict = {'Threshold': Thre,
              'Voxel Count': Vox_c1d,
              'Subject':Subjs,
              'Volume':Volume,
              'Data_type':DataType} 

pmaps_df = pd.DataFrame(data=pmaps_dict)

ary_vol = np.asarray(Volume).reshape(5,5)
#%%
ZaborPth = '/media/g/P01/Data/SPM_Anatomy/Zaborszky/4.2/'

os.chdir(ZaborPth)

Jubrain_v42=glob.glob('Ch-4_bl*')
Thre_Juv42 = []
Voxel_count_Juv42 = []
Vox_c1d_Juv42 = []
Atlas_Juv42 = []
for fileid in range(0,len(Jubrain_v42)):
     
     nii, hdrnii, affnii = fncLoadNii(Jubrain_v42[fileid])
     Num = []
     for th in threshold:
          num = np.count_nonzero(np.where(nii>th,upper,lower))
          Num.extend([num])
          Vox_c1d_Juv42.extend([num])
          Thre_Juv42.extend([th])
          if 'colin' in Jubrain_v42[fileid]:
               Atlas_Juv42.extend(['Jubrain MNI-Colin27 V4.2'])
          else:
               Atlas_Juv42.extend(['JuBrain MNI-ICBM152 V4.2'])
               
          
     Voxel_count_Juv42.append(Num)
     
Juv42_dict = {'Threshold': Thre_Juv42,
              'Volume': Vox_c1d_Juv42,
              'MNI Space':Atlas_Juv42} 

Juv42_df = pd.DataFrame(data=Juv42_dict)
Juv42_df_ICBM = Juv42_df.loc[Juv42_df['MNI Space']=='JuBrain MNI-ICBM152 V4.2']
Juv42_df_Cln  = Juv42_df.loc[Juv42_df['MNI Space']=='Jubrain MNI-Colin27 V4.2']
#%%
# JubrainV22: use a maximum probability map to delineate ch4
# ICBM[324], Colin27[279]
JubrainV22 = '/media/g/P01/Data/SPM_Anatomy/V22_Jubrain/MPM/'

os.chdir(JubrainV22)

Jubrain_v22=glob.glob('V22_bl*')
figPth = '/media/g/P01/Review/Figures/'




Voxel_count_Juv22 = []
Vox_c1d_Juv22 = []
Atlas_Juv22 = []
for fileid in range(0,len(Jubrain_v42)):
     
     nii, hdrnii, affnii = fncLoadNii(Jubrain_v22[fileid])
     num = np.count_nonzero(nii)
     
     Vox_c1d_Juv22.extend([num])
     if 'colin' in Jubrain_v42[fileid]:
          Atlas_Juv22.extend(['MNI-Colin27'])
     else:
          Atlas_Juv22.extend(['MNI-ICBM152'])
               
     
Juv22_dict = {
              'Volume': Vox_c1d_Juv22,
              'MNI Space':Atlas_Juv22} 

Juv22_df = pd.DataFrame(data=Juv22_dict)



#%% load zaborszky et al 2008 supplementary
excel=pd.read_excel('/media/g/P01/Data/SPM_Anatomy/Zaborszky/2008/Sheet2.xlsx')
excel['overlap'] = excel['overlap'].values / 100
excel['Ch4'] = excel['Ch4_l'].values + excel['ave_Ch4p'].values
excel_abpt3 = excel.loc[(excel['overlap']>0.2)]
excel_pt37 = excel_abpt3.loc[excel_abpt3['overlap']<0.8]
excel_pt37['Sources'] = ['Zaborszky et al., 2008' for i in range(len(excel_pt37))]
za2008 = excel_pt37[['overlap','Ch4','Sources']]
za2008_01=za2008.rename(columns={'Sources':'MNI Space','overlap':'Threshold','Ch4':'Volume'}) 
#%%
# Load nbM ROIs acquired from realiging 100um to 0.8 FatNav native space  
nbmROIPth = '/media/g/P01/Data/Bigbrain/nbmFatNav_seg/'
figPth = '/media/g/P01/Review/Figures/'
os.chdir(nbmROIPth)

nbmFatNavBB = glob.glob('*_bilateral_nbm.nii.gz')
nbmFatNavBB_left = glob.glob('*_left_nbm.nii.gz')
nbmFatNavBB_right = glob.glob('*_right_nbm.nii.gz')
Voxel_nbm_Fat = [] 
Volume_nbm_fav= []

for fileid in range(0,len(AllFiles)):
     
     nii, hdrnii, affnii = fncLoadNii(nbmFatNavBB[fileid])
     nii_l, hdrnii_l, affnii_l = fncLoadNii(nbmFatNavBB_left[fileid])
     nii_r, hdrnii_r, affnii_r = fncLoadNii(nbmFatNavBB_right[fileid])

          
     Voxel_nbm_Fat.append([ np.count_nonzero(nii), np.count_nonzero(nii_r),np.count_nonzero(nii_l)])
     Volume_nbm_fav.append(np.multiply([ np.count_nonzero(nii), np.count_nonzero(nii_r),np.count_nonzero(nii_l)],(0.8**3)))


ary_Volume_nbm_fav = np.asarray(Volume_nbm_fav)
ary_Vox_nbm_fav = np.asarray(Voxel_nbm_Fat)
ary_vol_fav_clm = ary_Volume_nbm_fav.T.reshape(-1)
ary_vox_fav_clm = ary_Vox_nbm_fav.T.reshape(-1)

Hemis = [i for i in ['BH','RH','LH'] for j in range(5)]
subjects = [j for i in range(3) for j in subjs_nme]
df_nbm_fatnav_bb = pd.DataFrame(data={'Volume':ary_vol_fav_clm,
                                      'Voxel':ary_vox_fav_clm,
                                      'Hemisphere':Hemis,
                                      'Subjects':subjects})



#%% Plot the results


lnewdth = 2.0
mksize = 5.0

mksize_scatter = 20

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

rc={'font.size': 9.0,
 'axes.labelsize': 11.0,
 'axes.titlesize': 12.0,
 'xtick.labelsize': 9.0,
 'ytick.labelsize': 9.0,
 'legend.fontsize': 9.0,
 'axes.linewidth': 1.25,
 'grid.linewidth': 1.0,
 'lines.linewidth': 2.0,
 'lines.markersize': 5.0,
 'patch.linewidth': 1.0,
 'xtick.major.width': 2,
 'ytick.major.width': 2,
 'xtick.minor.width': 1.0,
 'ytick.minor.width': 1.0,
 'xtick.major.size': 2.0,
 'ytick.major.size': 2.0,
 'xtick.minor.size': 1.0,
 'ytick.minor.size': 1.0,
 'legend.title_fontsize': 10.0}


left = 0.115  # the left side of the subplots of the figure
right = 0.983   # the right side of the subplots of the figure
bottom = 0.1  # the bottom of the subplots of the figure
top = 0.914     # the top of the subplots of the figure
wspace = 0.134  # the amount of width reserved for space between subplots,
              # expressed as a fraction of the average axis width
hspace = 0.2  # the amount of height reserved for space between subplots,


hisvol_halliday=[154,
134,
156,
128,
76]


hisvol_grinberg=99.06

#%% barplot of the nbm volume in individual space
# fig, ax0 = plt.subplots(1,1,figsize = [8,7], sharex = False,dpi=300)

# # plt.subplots_adjust(left = 0.115, right = 0.983,top = 0.914, wspace = 0.134 )
# plt.subplots_adjust(left = 0.07, right = 0.983,top = 0.914, wspace = 0.134 )
# plt.rcParams.update(**rc)

# sns.barplot(x='Subjects',y='Volume',hue='Hemisphere',ax=ax0,
#               data=df_nbm_fatnav_bb,palette='mako')

# ax0.set_ylabel("Volume ( " + 'mm\N{SUPERSCRIPT Three} )')
# plt.savefig(figPth + 'review_volume_report_fatnav_bigbrain_seg.pdf',dpi=300)
# plt.savefig(figPth + 'review_volume_report_fatnav_bigbrain_seg.png',dpi=300)
# plt.savefig(figPth + 'review_volume_report_fatnav_bigbrain_seg.svg',dpi=300)

#%%
fig, ax0 = plt.subplots(1,1,figsize = [8,7], sharex = False,dpi=300)

# plt.subplots_adjust(left = 0.115, right = 0.983,top = 0.914, wspace = 0.134 )
plt.subplots_adjust(left = 0.07, right = 0.983,top = 0.914, wspace = 0.134 )
plt.rcParams.update(**rc)



plt.plot(Juv42_df_ICBM['Threshold'], Juv42_df_ICBM['Voxel Count'], 'o:', color='grey',markersize=mksize, markerfacecolor='none',markeredgecolor='grey',markeredgewidth=1)
plt.plot(Juv42_df_Cln['Threshold'], Juv42_df_Cln['Voxel Count'], 'o--',color='grey',markersize=mksize, markerfacecolor='none',markeredgecolor='grey',markeredgewidth=1)
# plt.plot(Juv42_df_Cln['Threshold'], Juv42_df_Cln['Voxel Count'], 'o--',color='grey',markersize=mksize, markerfacecolor='none',markeredgecolor='grey',markeredgewidth=1)
plt.plot(excel_pt37['overlap'],excel_pt37['Ch4'], 'o-.', color='grey',markersize=mksize, markerfacecolor='none',markeredgecolor='grey',markeredgewidth=1)
std = np.std(ary_vol,axis=0)
plt.plot(np.asarray(threshold).reshape([-1]),np.mean(ary_vol,axis=0), 'o-', color='grey',markersize=mksize)
ax0.errorbar(np.asarray(threshold).reshape([-1]),np.mean(ary_vol,axis=0),yerr=std,color='grey',markersize=mksize, ecolor='grey',
       capsize=6)


pmaps_df_crop = pmaps_df.drop(['Voxel Count','Subject'],axis=1)


pmaps_df_crop_new = pmaps_df_crop.rename(columns={'Data_type':'MNI Space'}) 

cmb_df = pd.concat([Juv42_df,pmaps_df_crop_new,za2008_01],axis=0)

tmp_df = pd.DataFrame([[0.75, 154, 'Halliday et al., 1993'],
                       [0.75, 134, 'Halliday et al., 1993'],
                       [0.75, 156, 'Halliday et al., 1993'],
                       [0.75, 128, 'Halliday et al., 1993'],
                       [0.75, 76, 'Halliday et al., 1993'],
                       [0.75, 99.06, 'Gringberg et al., 2007'],
                       [0.8, 538, 'SPM Anatomy Toolbox v3.0'],
                       [0.8, 733, 'SPM Anatomy Toolbox v2.2'],
                       [0.85, 569, 'Jubrain MNI-ICBM152 v2.2'],
                       [0.85, 544, 'Jubrain MNI-Colin27 v2.2']], columns=['Threshold','Volume','MNI Space'])


cmb_df = pd.concat([Juv42_df,pmaps_df_crop_new,za2008_01],axis=0)




# sns.lineplot(x='Threshold',y='Volume',hue='Data_type',marker='o',markersize=mksize,ax=ax0,
#              data=pmaps_df,palette='cubehelix',linewidth = lnewdth,ci='sd',
#              err_style="bars",err_kws={'capsize':6})
# sns.lineplot(x='overlap',y='Ch4',hue='Sources',marker='o',markersize=mksize,ax=ax0,
#              data=excel_pt37,palette='rocket',linewidth = lnewdth,ci='sd')


hloc= 0.8
step=0.1


     

fig, ax0 = plt.subplots(1,1,figsize = [8,7],dpi=300)

plt.rcParams.update(**rc)
plt.subplots_adjust(left = 0.115, right = 0.983,top = 0.914, wspace = 0.134 )
for val in hisvol_halliday:
     # ax0.scatter(0.75,val,c='darksalmon', label='Halliday et al., 1993')
     if val == 154:
          
          ax0.scatter(hloc,val,s=mksize_scatter,edgecolors='black', facecolors='none',label='Halliday et al., 1993')
     else:
          ax0.scatter(hloc,val,s=mksize_scatter,edgecolors='black', facecolors='none')
ax0.scatter(hloc,hisvol_grinberg,c='black',s=mksize_scatter,label='Gringberg et al., 2007')          
plt.plot(Juv42_df_ICBM['Threshold'], Juv42_df_ICBM['Volume'], 'o:', color='grey',marker=">",markersize=mksize,markeredgecolor='grey',markeredgewidth=1,label='MNI-ICBM152 v4.2')
plt.plot(Juv42_df_Cln['Threshold'], Juv42_df_Cln['Volume'], 'o--',color='grey',marker="P", markersize=mksize,markeredgecolor='grey',markeredgewidth=1,label='MNI-Colin27 v4.2')
# plt.plot(Juv42_df_Cln['Threshold'], Juv42_df_Cln['Voxel Count'], 'o--',color='grey',markersize=mksize, markerfacecolor='none',markeredgecolor='grey',markeredgewidth=1)
plt.plot(excel_pt37['overlap'],excel_pt37['Ch4'], 'o-.', color='grey',marker="D", markersize=mksize,markeredgecolor='grey',markeredgewidth=1,label='Zaborszky et al., 2008')
std = np.std(ary_vol,axis=0)
plt.plot(np.asarray(threshold).reshape([-1]),np.mean(ary_vol,axis=0), 'o-', color='grey',markersize=mksize,label='in vivo native space')
ax0.errorbar(np.asarray(threshold).reshape([-1]),np.mean(ary_vol,axis=0),yerr=std,color='grey',markersize=mksize, ecolor='grey',
       capsize=6)

ax0.scatter(hloc+step,538,c='darkgray',s=mksize_scatter, label='v3.0')
ax0.scatter(hloc+step,733,edgecolors='darkgray',s=mksize_scatter, facecolors='none',label='v2.2')
ax0.scatter(hloc+step*2,569,c='slategray',s=mksize_scatter, label='Jubrain MNI-ICBM152 v2.2')
ax0.scatter(hloc+step*2,544,s=mksize_scatter, edgecolors='slategray', facecolors='none',label='Jubrain MNI-Colin27 v2.2')
ax0.scatter(0.4,1377,s=mksize_scatter,c='olive')
ax0.annotate('Markello et al., 2018', (0.403, 1377),fontsize=9.0)

ax0.annotate('Histology', (0.503, 1000),fontsize=10.0)
ax0.annotate('Probabilistic Maps', (0.503, 900),fontsize=10.0)
ax0.annotate('Maximum Probability Maps', (0.503, 800),fontsize=10.0)
ax0.annotate('Jubrain', (0.303, 800),fontsize=9.0)
# ax0.annotate('Jubrain v2.2', (0.603, 900),fontsize=10.0)
# ax0.annotate('Jubrain v4.2', (0.403, 800),fontsize=10.0)
handles,_ = ax0.get_legend_handles_labels()
# labels = ['LH_BigBrain_NBM','RH_BigBrain_NBM','LH_ICBM_NBM','LH_Colin_NBM','RH_ICBM_NBM','RH_Colin_NBM']
# labels = ['Jubrain MNI-ICBM152 v4.2','Jubrain MNI-Colin27 v4.2','in vivo 7T 0.8 mm iso',
#           'Gringberg et al., 2007','Halliday et al., 1993','Halliday et al., 1993',
#           'Halliday et al., 1993','Halliday et al., 1993','Halliday et al., 1993',
#           'Zaborszky et al., 2008','Jubrain MNI-ICBM152 v2.2','Jubrain MNI-Colin27 v2.2',
#           'SPM Anatomy Toolbox v3.0','SPM Anatomy Toolbox v2.2']

# labels = ['Jubrain MNI-ICBM152 v4.2','Jubrain MNI-Colin27 v4.2','in vivo 7T 0.8 mm iso',
#            'Zaborszky et al., 2008','Gringberg et al., 2007','Halliday et al., 1993',
#          'Jubrain MNI-ICBM152 v2.2','Jubrain MNI-Colin27 v2.2',
#           'SPM Anatomy Toolbox v3.0','SPM Anatomy Toolbox v2.2']

labels = ['MNI-ICBM152 v4.2','MNI-Colin27 v4.2','in vivo 7T 0.8 mm iso',
           'Zaborszky et al., 2008','Halliday et al., 1993','Gringberg et al., 2007',
         'SPM Anat. Tbx v3.0','SPM Anat. Tbx v2.2','MNI-ICBM152 v2.2','MNI-Colin27 v2.2',
          ]


ax0.legend(handles=handles[0:],labels = labels)
# plt.legend(loc=1)
# plt.gca().legend(loc='center left', bbox_to_anchor=(1.0, 0.5),borderaxespad=0.)
xticklabels = np.round(ax0.get_xticks(),2).astype(str)

xticklabels[-3:] = ['Histology','SPM Anat. Tbx','JuBrain V2.2']
# xlables = np.hstack([xticklabels,])
ax0.set_xticklabels(xticklabels[1:],rotation=15)

ax0.set_xlabel("Threshold of probabilistic maps")
ax0.set_ylabel("Volume ( " + 'mm\N{SUPERSCRIPT Three} )')
fig.suptitle("nbM Volume from Published Literature and In vivo 7T Datasets")
plt.savefig(figPth + 'review_volume_report_lt_invivo_01.pdf',dpi=300)
plt.savefig(figPth + 'review_volume_report_lt_invivo_01.png',dpi=300)
plt.savefig(figPth + 'review_volume_report_lt_invivo_01.svg',dpi=300)







hloc= 0.8
step=0.1
fig, ax0 = plt.subplots(1,1,figsize = [8,7], sharex = False,dpi=300)
plt.rcParams.update(**rc)
plt.subplots_adjust(left = 0.115, right = 0.983,top = 0.914, wspace = 0.134 )
for val in hisvol_halliday:
     # ax0.scatter(0.75,val,c='darksalmon', label='Halliday et al., 1993')
     if val == 154:
          
          ax0.scatter(hloc,val,s=mksize_scatter,edgecolors='black', facecolors='none',label='Halliday et al., 1993')
     else:
          ax0.scatter(hloc,val,s=mksize_scatter,edgecolors='black', facecolors='none')
ax0.scatter(hloc,hisvol_grinberg,c='black',s=mksize_scatter,label='Gringberg et al., 2007')          
plt.plot(Juv42_df_ICBM['Threshold'], Juv42_df_ICBM['Volume'],  'o:', color='grey',marker=">", markersize=mksize,markeredgecolor='grey',markeredgewidth=1,label='MNI-ICBM152 v4.2')
plt.plot(Juv42_df_Cln['Threshold'], Juv42_df_Cln['Volume'],'o--', color='grey',marker="H", markersize=mksize,markeredgecolor='grey',markeredgewidth=1,label='MNI-Colin27 v4.2')
# plt.plot(Juv42_df_Cln['Threshold'], Juv42_df_Cln['Voxel Count'], 'o--',color='grey',markersize=mksize, markerfacecolor='none',markeredgecolor='grey',markeredgewidth=1)
plt.plot(excel_pt37['overlap'],excel_pt37['Ch4'], 'o-.', color='grey',marker="D", markersize=mksize,markeredgecolor='grey',markeredgewidth=1,label='Zaborszky et al., 2008')
std = np.std(ary_vol,axis=0)
plt.plot(np.asarray(threshold).reshape([-1]),np.mean(ary_vol,axis=0), 'o-', color='grey', marker="s", markersize=mksize,label='in vivo native space')
ax0.errorbar(np.asarray(threshold).reshape([-1]),np.mean(ary_vol,axis=0),yerr=std,color='grey',markersize=mksize, ecolor='grey',
       capsize=6)

ax0.scatter(hloc+step,538,c='darkgrey',s=mksize_scatter, label='SPM Anatomy Toolbox v3.0')
ax0.scatter(hloc+step,733,edgecolors='darkgrey',s=mksize_scatter, facecolors='none',label='SPM Anatomy Toolbox v2.2')
ax0.scatter(hloc+step*2,569,c='darkgrey',s=mksize_scatter, label='Jubrain MNI-ICBM152 v2.2')
ax0.scatter(hloc+step*2,544,s=mksize_scatter, edgecolors='darkgrey', facecolors='none',label='Jubrain MNI-Colin27 v2.2')
ax0.scatter(0.4,1377,s=mksize_scatter,c='darkgrey')
ax0.annotate('Markello et al., 2018', (0.404, 1370),fontsize=9.0)

ax0.scatter(0.5,791,s=mksize_scatter,c='darkgrey')
ax0.annotate('Yuan et al., 2019', (0.504, 780),fontsize=9.0)


ax0.scatter(0.5,240,s=mksize_scatter,c='darkgrey')
ax0.vlines(0.5,200,300,colors='darkgrey',linestyles='dotted')
# ax0.hlines(200,0.49,0.51,colors='k')
ax0.hlines(300,0.495,0.505,colors='darkgrey')
ax0.hlines(200,0.495,0.505,colors='darkgrey')
ax0.annotate('Fernández-Cabello et al. 2020', (0.504, 230),fontsize=9.0)

# ax0.scatter(0.5,250,s=50,facecolors='none',edgecolors='dimgray')

# ax0.annotate('Jubrain v2.2', (0.683, 900),fontsize=10.0)
# ax0.annotate('Jubrain v4.2', (0.403, 800),fontsize=10.0)
handles,_ = ax0.get_legend_handles_labels()

labels = ['Jubrain MNI-ICBM152 v4.2','Jubrain MNI-Colin27 v4.2','in vivo 7T 0.8 mm iso',
           'Zaborszky et al., 2008','Gringberg et al., 2007','Halliday et al., 1993',
         'Jubrain MNI-ICBM152 v2.2','Jubrain MNI-Colin27 v2.2',
          'SPM Anatomy Toolbox v3.0','SPM Anatomy Toolbox v2.2']


ax0.legend(handles=handles[0:],labels = labels)
plt.legend(loc=1)

xticklabels = np.round(ax0.get_xticks(),2).astype(str)

xticklabels[-4:-1] = ['Histology','SPM Anat. Tbx','JuBrain V2.2']
ax0.set_xticklabels(xticklabels[0:],rotation=15)

ax0.set_xlabel("Threshold of probabilistic maps")
ax0.set_ylabel("Volume ( " + 'mm\N{SUPERSCRIPT Three} )')
# fig.suptitle("nbM Volume from Published Literature and In vivo 7T Datasets")

plt.savefig(figPth + 'review_volume_report_lt_00.pdf',dpi=300)
plt.savefig(figPth + 'review_volume_report_lt_00.png',dpi=300)
plt.savefig(figPth + 'review_volume_report_lt_00.svg',dpi=300)


