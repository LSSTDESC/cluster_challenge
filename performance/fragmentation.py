#!/usr/bin/env python
# coding: utf-8

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from astropy.table import Table
import sys
import os

###clevar
import clevar
from clevar.catalog import ClCatalog

#
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
else :
    print('careful, wrong option', str(sys.argv[2]), 'not allowed')

outpath += matching + '_matching/'     
    
c1 = ClCatalog.read(cat_path + 'c1_' + matching + '.fits', 'c1', full = True)
c2 = ClCatalog.read(cat_path + 'c2_' + matching + '.fits', 'c2', full = True)
c_merged_12 = ClCatalog.read(cat_path+'output_catalog_' + matching + '.fits', 'merged',  full=True)
c_merged_12 = c_merged_12[[type(siu) != ma.core.MaskedConstant for siu in c_merged_12['cat2_mt_multi_other']]] #Prevent from mask errors
# restrict to z<1.15
# c1=c1[c1.data['z']<1.15]

# restrict to matched pairs
# mt1, mt2 = get_matched_pairs(c1, c2, 'cross', None, None) 

# plot style
figx=10
figy=7

print('------Fragmentation------')
nbins_x = 9
bin1 = np.linspace(13, 15, nbins_x)
bin2 = [.2,.5,.8,1,1.2,1.8]
labels=['0.2-0.5','0.5-0.8','0.8-1.0','1.0-1.2','1.2-1.8']
plt.xlim(13,15)
nbins_x = nbins_x-1
bin_range = [min(bin1), max(bin1)]
frag = np.empty([len(bin2),nbins_x])
bin_x = np.empty([nbins_x])

for ix in range(nbins_x):
     bin_x[ix] = 0.5 * (bin1[ix] + bin1[ix+1])

for i in range(0,len(bin2)-1):
    cut1 = bin2[i]
    cut2 = bin2[i+1]
    filter1 = np.logical_and(c_merged_12.data['cat1_z'] > cut1, c_merged_12.data['cat1_z'] < cut2) #Attention comp avc z_cl
    c_clusters_matched = c_merged_12[filter1]
    #print(c_clusters_matched)
    filter2 = np.logical_and(c2.data['z'] > cut1, c2.data['z'] < cut2)
    c_clusters = c2.data[filter2]
    #print(c_clusters)
    h_r_clusters_matched = np.histogram(c_clusters_matched[[len(siu.split(','))>=2 for siu in c_clusters_matched['cat2_mt_multi_other']]]['cat2_log_mass'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    h_r_clusters  = np.histogram(c_clusters['log_mass'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    #print(h_r_clusters_matched)
    #print(h_r_clusters)
    frag[i] = np.divide(h_r_clusters_matched[0],h_r_clusters[0],where=(h_r_clusters[0]!=0))
    for j in range(len(frag[i])):
        if h_r_clusters_matched[0][j]<10 or h_r_clusters[0][j]<10:
            frag[i][j] = np.nan
    plt.scatter(bin_x, frag[i], label=labels[i], marker= ".", s=30)
    plt.plot(bin_x, frag[i])


#plot in bins of redshift
plt.ylabel('Fragmentation')
plt.xlabel('log(M/$M_{\odot}$)', fontsize = 13)
plt.legend()
plt.ylim(10**-4,10)
plt.yscale('log')
plt.xlim(13.0,15.5) 
plt.title('Amico-cosmoDC2')
plt.savefig(outpath+'nobckgfragmentation.png', bbox_inches='tight', format='png', transparent=True)
