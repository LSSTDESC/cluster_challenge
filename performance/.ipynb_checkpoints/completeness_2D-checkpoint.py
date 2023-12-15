#!/usr/bin/env python
# coding: utf-8

###import
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import sys
import os

###clevar
import clevar
from clevar.catalog import ClCatalog

inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/'
cat_path = inpath + str(sys.argv[1])

print('Using catalogs located at :', cat_path)

# #---Select catalogs to match---##

outpath = '/pbs/home/n/namourou/test_jupyter/cluster_challenge/plots/'
outpath += str(sys.argv[1])

print('Saving plots at :', outpath)

if not os.path.exists(outpath):
    os.makedirs(outpath)

##########select case
if str(sys.argv[2]) == 'p_matching':
    matching = 'p'
elif str(sys.argv[2]) == 'mb_matching':
    matching = 'mb'

outpath += matching + '_matching/'    
    
c1 = ClCatalog.read(cat_path + 'c1_' + matching + '.fits', 'c1', full = True)
c2 = ClCatalog.read(cat_path + 'c2_' + matching + '.fits', 'c2', full = True)
c_merged_12 = ClCatalog.read(cat_path+'output_catalog_' + matching + '.fits', 'merged',  full=True)


#plot style
figx=10
figy=7

#2D : z-M plane
#nbins_x = 18 #AM
nbins_x = 10 #RM
nbins_y = 21
x_bins = np.linspace(0,1.8,nbins_x)
y_bins = np.logspace(13,15,nbins_y)
xbin_range = [min(x_bins), max(x_bins)]
ybin_range = [min(y_bins), max(y_bins)]
h2_z_halos_matched = np.histogram2d(c_merged_12['cat2_z'], c_merged_12['cat2_mass'], bins=(x_bins, y_bins), 
                                    range=[xbin_range, ybin_range], normed=None, weights=None, density=None)
h2_z_halos = np.histogram2d(c2.data['z'], c2.data['mass'], bins=(x_bins, y_bins),
                            range=[xbin_range, ybin_range], normed=None, weights=None, density=None)
number_of_match = h2_z_halos_matched[0]
number_of_halo = h2_z_halos[0]
compl_2d = np.divide(number_of_match, number_of_halo, where=(number_of_halo!=0))
fig, ax = plt.subplots(figsize =(12,6))
x, y = np.meshgrid(x_bins, y_bins)
#to avoid seeing 0 values
compl_2d[compl_2d==0] = np.nan
#print(x)
#print(y)
c = ax.pcolormesh(x, y, compl_2d.T, cmap='jet', vmin=0, vmax=1)
ax.set_xlim(0,2.0)
ax.set_ylim(10**13,10**15)
ax.set_xlabel('z', fontsize = 13)
ax.set_ylabel('halo_mass', fontsize = 13)
ax.set_yscale('log')
ax.set_title('2D plot completness', fontsize = 13)
fig.colorbar(c, ax=ax, label = 'Completness')
plt.savefig(outpath + 'compl2D_plot_' +  matching + '.png')


