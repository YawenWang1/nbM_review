#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 17:44:45 2021
This script is to plot histogram after the simulation on Bigbrain dataset

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

Pth = '/media/g/P01/Yawen/tmp_msk_files/'

FigPth = '/media/g/P01/Review/Figures/'

df = pd.read_excel(Pth + 'nbm_bigbrain_simul_ds.xlsx')

#%% Plot all the hisgram in one figure
my_dpi = 125
rc_hist={'font.size': 12.0,
 'axes.labelsize': 12.0,
 'axes.titlesize': 12.0,
 'xtick.labelsize': 12.0,
 'ytick.labelsize': 12.0,
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

# plot bilateral nbM volume
## Separate data for left and right hemispheres
df_l = df.loc[df['Hemisphere'] == 'Left']   

df_r = df.loc[df['Hemisphere'] == 'Right']   

## Caluclate bilateral data
df_bl = df_l.copy()

df_bl[df.columns[2]] = 'Bilteral'
df_bl[df.columns[3]] = df_l[df.columns[3]].values + df_r[df.columns[3]].values
df_bl[df.columns[4]] = df_l[df.columns[4]].values + df_r[df.columns[4]].values
boundaries = np.linspace(0,150,10) 
fig, axs = plt.subplots(6,5,figsize=[11,8],sharex = True,sharey=True,dpi=my_dpi)

plt.rcParams.update(**rc_hist)
plt.ylim([0, 5000])
plt.xlim([0, 150])
plt.subplots_adjust(left =0.05, right = 0.985,top = 0.975, wspace = 0.134, bottom=0.03 )
for rowid in range(6):
     
     for clnid in range(5):
          
          curres = rowid * 5  + clnid  
     
          axs[rowid, clnid].hist(df_bl.loc[df['Resolution']==res[curres]]['Volume Size'], bins=boundaries, label='histogram')
          # axs[rowid, clnid].set_title('R : {}  mm iso.'.format(res[curres]))
          
          axs[rowid, clnid].annotate('R:{}  mm '.format(res[curres]), (80, 3500),fontsize=12.0)
     

plt.savefig(FigPth + 'hist res ' + str(0) + ' ' + str(30) + '_bl_v3_bins_10.pdf',dpi=300)
plt.savefig(FigPth + 'hist res ' + str(0) + ' ' + str(30) + '_bl_v3_bins_10.svg',dpi=300)
plt.savefig(FigPth + 'hist res ' + str(0) + ' ' + str(30) + '_bl_v3_bins_10.png',dpi=300)
plt.close('all')








     