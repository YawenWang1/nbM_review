#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 17:44:45 2021
This script is to plot histogram after the simulation

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
from scipy.optimize import curve_fit
#%% Load relevant files 
Pth = '/media/g/P01/Yawen/tmp_msk_files/'

FigPth = '/media/g/P01/Review/Figures/'

#-----------------------------------------------------------------------------
# Combine all the xlsx files
# df = pd.DataFrame()

# for df_id in df_files:
     
#      curr_df = pd.read_excel(df_id)
#      df = pd.concat([df,curr_df])
#      # print(len(curr_df), df_id)
     
# df = pd.concat([df,df_orig])
# df.to_excel(Pth + 'nbm_bigbrain_simul_ds.xlsx')
#------------------------------------------------------------------------------

df = pd.read_excel(Pth + 'nbm_bigbrain_simul_ds.xlsx')
## Separate data for left and right hemispheres
df_l = df.loc[df['Hemisphere'] == 'Left']   

df_r = df.loc[df['Hemisphere'] == 'Right']   

## Caluclate bilateral data
df_bl = df_l.copy()

df_bl[df.columns[2]] = 'Bilteral'
df_bl[df.columns[3]] = df_l[df.columns[3]].values + df_r[df.columns[3]].values
df_bl[df.columns[4]] = df_l[df.columns[4]].values + df_r[df.columns[4]].values

#%% Calculate mean and coefficiant of variant for each resolution
res = df['Resolution'].unique()

mean_bl, mean_r, mean_l, cv_bl, cv_r, cv_l  = [], [], [], [], [], []

for resid in range(len(res)):
     curr_df_bl = df_bl.loc[df_bl['Resolution']==res[resid]]['Volume Size']
     
     mean_bl.append(curr_df_bl.mean())
     
     cv_bl.append(np.divide(curr_df_bl.std(),curr_df_bl.mean()))
     
     
     curr_df_r = df_r.loc[df_r['Resolution']==res[resid]]['Volume Size']
     
     mean_r.append(curr_df_r.mean())
     
     cv_r.append(np.divide(curr_df_r.std(),curr_df_r.mean()))
     
     
     curr_df_l = df_l.loc[df_l['Resolution']==res[resid]]['Volume Size']
     
     mean_l.append(curr_df_l.mean())
     
     cv_l.append(np.divide(curr_df_l.std(),curr_df_l.mean()))
     
# Create a dataframe      
CV    = np.asarray(cv_bl + cv_r + cv_l)
Mean  = np.asarray(mean_bl + mean_r + mean_l)
Res   = np.asarray([j for i in range(3) for j in res])
Hemis = np.asarray([i for i in ['BH', 'RH', 'LH'] for j in range(len(res))])


CV    = cv_bl + cv_r + cv_l
Mean  = mean_bl + mean_r + mean_l
Res   = [j for i in range(3) for j in res]
Hemis = [i for i in ['BH', 'RH', 'LH'] for j in range(len(res))]

df_mcv = pd.DataFrame(data={'Coefficient of variation of volume': CV,
                            'Volume mean': Mean,
                            'Resolution': Res,
                            'Hemisphere':Hemis})

df_mcv.to_excel(FigPth + 'Mean_CV_pt_01_30.xlsx')
#%% Mean and CV across resolution segmentations
df_partial_submm = df_mcv.loc[df['Resolution'] < 1.1]
tmp_bh = df_mcv[((df_mcv['Resolution'] < 1.1) & (df_mcv['Hemisphere'] == 'BH')) ]
tmp_lh = df_mcv[((df_mcv['Resolution'] < 1.1) & (df_mcv['Hemisphere'] == 'LH')) ]
tmp_rh = df_mcv[((df_mcv['Resolution'] < 1.1) & (df_mcv['Hemisphere'] == 'RH')) ]

print(tmp_bh['Coefficient of variation of volume'].mean())
print(tmp_bh['Volume mean'].mean())
print(tmp_bh['Coefficient of variation of volume'].std())
print(tmp_bh['Volume mean'].std())
print(tmp_lh['Coefficient of variation of volume'].mean())
print(tmp_lh['Volume mean'].mean())
print(tmp_lh['Coefficient of variation of volume'].std())
print(tmp_lh['Volume mean'].std())
print(tmp_rh['Coefficient of variation of volume'].mean())
print(tmp_rh['Volume mean'].mean())
print(tmp_rh['Coefficient of variation of volume'].std())
print(tmp_rh['Volume mean'].std())


# 0.14692438745353698
# 59.89620094984167
# 0.00816168361879775
# 1.4305585758023214

# 0.1033097789414212
# 29.463081286452223
# 0.03581721458331479
# 0.9096146408495623

# 0.2218988420660743
# 30.437185652391314
# 0.008450609568404675
# 0.5609580652180576

df_partial_abvmm_01 = df_mcv.loc[df['Resolution'] > 1.0]  

df_partial_abv12 = df_partial_abvmm_01.loc[df_partial_abvmm_01['Resolution'] < 2.1] 
tmp_bh_12 = df_partial_abv12[((df_partial_abv12['Resolution'] > 1.0) & (df_partial_abv12['Hemisphere'] == 'BH')) ]
tmp_lh_12 = df_partial_abv12[((df_partial_abv12['Resolution'] > 1.0) & (df_partial_abv12['Hemisphere'] == 'LH')) ]
tmp_rh_12 = df_partial_abv12[((df_partial_abv12['Resolution'] > 1.0) & (df_partial_abv12['Hemisphere'] == 'RH')) ]



print(tmp_bh_12['Coefficient of variation of volume'].mean())
print(tmp_bh_12['Volume mean'].mean())
print(tmp_bh_12['Coefficient of variation of volume'].std())
print(tmp_bh_12['Volume mean'].std())
print('/n')
print(tmp_lh_12['Coefficient of variation of volume'].mean())
print(tmp_lh_12['Volume mean'].mean())
print(tmp_lh_12['Coefficient of variation of volume'].std())
print(tmp_lh_12['Volume mean'].std())
print('/n')
print(tmp_rh_12['Coefficient of variation of volume'].mean())
print(tmp_rh_12['Volume mean'].mean())
print(tmp_rh_12['Coefficient of variation of volume'].std())
print(tmp_rh_12['Volume mean'].std())


# 0.255803085990428
# 53.61139096817227
# 0.08844399811414587
# 4.279903913077699
# /n
# 0.3432434781775044
# 25.748558540243426
# 0.17742304356856803
# 2.916944427338156
# /n
# 0.36860075489946675
# 27.8651320446594
# 0.10446843338787024
# 2.2480427355933696

df_partial_abv23 = df_partial_abvmm_01.loc[df_partial_abvmm_01['Resolution'] > 2.0] 

tmp_bh_23 = df_partial_abv23[df_partial_abv23['Hemisphere'] == 'BH']
tmp_lh_23 = df_partial_abv23[df_partial_abv23['Hemisphere'] == 'LH']
tmp_rh_23 = df_partial_abv23[df_partial_abv23['Hemisphere'] == 'RH']


print(tmp_bh_23['Coefficient of variation of volume'].mean())
print(tmp_bh_23['Volume mean'].mean())
print(tmp_bh_23['Coefficient of variation of volume'].std())
print(tmp_bh_23['Volume mean'].std())
print('/n')
print(tmp_lh_23['Coefficient of variation of volume'].mean())
print(tmp_lh_23['Volume mean'].mean())
print(tmp_lh_23['Coefficient of variation of volume'].std())
print(tmp_lh_23['Volume mean'].std())
print('/n')
print(tmp_rh_23['Coefficient of variation of volume'].mean())
print(tmp_rh_23['Volume mean'].mean())
print(tmp_rh_23['Coefficient of variation of volume'].std())
print(tmp_rh_23['Volume mean'].std())

# 0.5804328324012018
# 50.34651329778283
# 0.1362360608705913
# 6.933103000868983
# /n
# 0.9727705180650329
# 23.43216452257924
# 0.30919020443473244
# 6.615960535406939
# /n
# 0.7362287982621599
# 26.907983169471464
# 0.11761479783399315
# 1.7607118271265474


tmp_bh_10 = df_bl[df_bl['Resolution'] == 1]
tmp_lh_10 = df_l[df_l['Resolution'] == 1]
tmp_rh_10 = df_r[df_r['Resolution'] == 1]
print('The mean volume at 1.0 mm iso is {} and std is {}'.format(tmp_bh_10['Num of Voxels'].mean(),tmp_bh_10['Num of Voxels'].std()))
# The mean volume at 1.0 mm iso is 57.085985669055155 and std is 9.392582458113504

tmp_bh_20 = df_bl[df_bl['Resolution'] == 2]
tmp_lh_20 = df_l[df_l['Resolution'] == 2]
tmp_rh_20 = df_r[df_r['Resolution'] == 2]
print('The mean volume at 1.0 mm iso is {} and std is {}'.format(tmp_bh_20['Num of Voxels'].mean(),tmp_bh_20['Num of Voxels'].std()))
# The mean volume at 1.0 mm iso is 5.9928345275787365 and std is 2.4266125870125297

tmp_bh_30 = df_bl[df_bl['Resolution'] == 3]
tmp_lh_30 = df_l[df_l['Resolution'] == 3]
tmp_rh_30 = df_r[df_r['Resolution'] == 3]
print('The mean volume at 1.0 mm iso is {} and std is {}'.format(tmp_bh_30['Num of Voxels'].mean(),tmp_bh_30['Num of Voxels'].std()))
# The mean volume at 1.0 mm iso is 1.698883519413431 and std is 1.2879494373861684
#%% set up for the plot
rc_hist={'font.size': 12.0,
 'axes.labelsize': 17.0,
 'axes.titlesize': 17.0,
 'xtick.labelsize': 17.0,
 'ytick.labelsize': 17.0,
 'legend.fontsize': 17.0,
 'axes.linewidth': 1.25,
 'grid.linewidth': 1.0,
 'lines.linewidth': 4,
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
 'legend.title_fontsize': 17.0,
 'font.family':'Arial'}


#%%

plt.figure(figsize=(11.7,8.3),dpi=300)
plt.rcParams.update(**rc_hist)
plt.subplots_adjust(left = 0.065, right = 0.983,top = 0.974, wspace = 0.134, bottom=0.08 )   
ax0 = sns.barplot(x='Resolution',y='Coefficient of variation of volume',hue='Hemisphere',data=df_mcv,palette='mako')
ax0.set_xlabel("Spatial resolution ( mm ) isotropic")       
plt.savefig(FigPth + 'CV pt' + str(1) + ' ' + str(30) + ' bar.pdf',dpi=300)
plt.savefig(FigPth + 'CV pt' + str(1) + ' ' + str(30) + ' bar.svg',dpi=300)
plt.savefig(FigPth + 'CV pt' + str(1) + ' ' + str(30) + ' bar.png',dpi=300)
plt.close('all')



# plt.figure(figsize=(11.7,8.3),dpi=300)
# plt.rcParams.update(**rc_hist)
# plt.subplots_adjust(left = 0.065, right = 0.983,top = 0.974, wspace = 0.134, bottom=0.08 )   
# ax1 = sns.lineplot(x='Resolution',y='Coefficient of variation of volume',hue='Hemisphere',data=df_mcv,palette='mako')
# ax1.set_xlabel("Spatial resolution ( mm ) isotropic")        
# plt.savefig(FigPth + 'CV pt' + str(1) + ' ' + str(30) + ' line.pdf',dpi=300)
# plt.savefig(FigPth + 'CV pt' + str(1) + ' ' + str(30) + ' line.svg',dpi=300)
# plt.savefig(FigPth + 'CV pt' + str(1) + ' ' + str(30) + ' line.png',dpi=300)
# plt.close('all')


plt.figure(figsize=(11.7,8.3),dpi=300)
plt.rcParams.update(**rc_hist)
plt.subplots_adjust(left = 0.065, right = 0.983,top = 0.974, wspace = 0.134, bottom=0.08 )      
ax2 = sns.barplot(x='Resolution',y='Volume mean',hue='Hemisphere',data=df_mcv,palette='mako')
ax2.set_xlabel("Spatial resolution ( mm ) isotropic")           
plt.savefig(FigPth + 'Mean pt' + str(1) + ' ' + str(30) + ' bar.pdf',dpi=300)
plt.savefig(FigPth + 'Mean pt' + str(1) + ' ' + str(30) + ' bar.svg',dpi=300)
plt.savefig(FigPth + 'Mean pt' + str(1) + ' ' + str(30) + ' bar.png',dpi=300)
plt.close('all')



