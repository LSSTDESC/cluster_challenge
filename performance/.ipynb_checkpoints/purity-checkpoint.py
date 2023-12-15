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

outpath += matching + '_matching/'     
    
c1 = ClCatalog.read(cat_path + 'c1_' + matching + '.fits', 'c1', full = True)
c2 = ClCatalog.read(cat_path + 'c2_' + matching + '.fits', 'c2', full = True)
c_merged_12 = ClCatalog.read(cat_path+'output_catalog_' + matching + '.fits', 'merged',  full=True)


#plot style

figx=10
figy=7
print('Purity vs mass in redshift bins')
bin_range = [0,100]
nbins_x = 10
zbins = [0,0.5,0.8,1.0,1.2,1.8] #Amico
purity_z_raw = np.empty([len(zbins),nbins_x])

bin_x = np.empty([nbins_x])
x_bins = np.linspace(0,100,nbins_x+1)
labels=['0-0.5','0.5-0.8','0.8-1.0','1.0-1.2', '1.2-1.8']

for ix in range(nbins_x):
     bin_x[ix] = 0.5 * (x_bins[ix] + x_bins[ix+1])
        
for i in range(0,len(zbins)-1):
    cut1 = zbins[i]
    cut2 = zbins[i+1]
    filter1 = np.logical_and(c_merged_12.data['cat1_z'] > cut1, c_merged_12.data['cat1_z'] < cut2)
    c_clusters_matched = c_merged_12[filter1]
    #print(c_clusters_matched)
    filter2 = np.logical_and(c1.data['z'] > cut1, c1.data['z'] < cut2)
    c_clusters = c1.data[filter2]
    #print(c_clusters)
    h_r_clusters_matched = np.histogram(c_clusters_matched['cat1_mass'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    h_r_clusters  = np.histogram(c_clusters['mass'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    #print(h_r_clusters_matched)
    #print(h_r_clusters)
    purity_z_raw[i] = np.divide(h_r_clusters_matched[0],h_r_clusters[0],where=(h_r_clusters[0]!=0))
    for j in range(len(purity_z_raw[i])):
        if h_r_clusters_matched[0][j]<10 or h_r_clusters[0][j]<10:
            purity_z_raw[i][j] = np.nan
    plt.scatter(bin_x, purity_z_raw[i], label=labels[i], marker= ".", s=30)
    plt.plot(bin_x, purity_z_raw[i])
#plot in bins of redshift
plt.xlabel('$\lambda^*$', fontsize = 13)
plt.ylabel('Purity', fontsize = 13)
plt.legend()
plt.ylim(0,1.2)
plt.xlim(0,110) 
plt.title('AMICO-cosmoDC2', fontsize = 13)
plt.savefig(outpath + 'nobckgpurity_' + matching + '.png', bbox_inches='tight', format='png', transparent=True)