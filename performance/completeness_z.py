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

##---Select catalogs to match---##

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
c_merged_12 = ClCatalog.read(cat_path + 'output_catalog_' + matching + '.fits', 'merged',  full=True)


#plot style
figx=10
figy=7

print('completeness vs z in mass bins') 
param1, param2, param3, param4 = 'cat2_z', 'z', 'cat2_mass', 'mass'
nbins_x = 10
bin1 = np.linspace(0.0, 1.8, nbins_x) #For AMICO
bin2 = [10**13,10**13.5,10**14,10**14.3,10**15]
labels=['$10^{13}-10^{13.5}$','$10^{13.5}-10^{14}$', '$10^{14}-10^{14.3}$','$10^{14.3}-10^{15}$']
plt.xlim(0,2.0) #For AMICO

nbins_x -= 1
bin_range = [min(bin1), max(bin1)]
compl = np.empty([len(bin2),nbins_x])
bin_x = np.empty([nbins_x])
for ix in range(nbins_x):
    bin_x[ix] = 0.5 * (bin1[ix] + bin1[ix+1])
for i in range(0,len(bin2)-1):
    cut1 = bin2[i]
    cut2 = bin2[i+1]
    filter1 = np.logical_and(c_merged_12.data[param3] > cut1, c_merged_12.data[param3] < cut2)
    filter2 = np.logical_and(c2.data[param4] > cut1, c2.data[param4] < cut2)
    c_halos_matched = c_merged_12[filter1]
    c_halos = c2.data[filter2]
    h_r_halos_matched = np.histogram(c_halos_matched[param1], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    h_r_halos  = np.histogram(c_halos[param2], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    compl[i] = np.divide(h_r_halos_matched[0], h_r_halos[0], where=(h_r_halos[0]!=0))
    for j in range(len(compl[i])):
        if h_r_halos_matched[0][j]<10 or h_r_halos[0][j]<10:
            compl[i][j] = np.nan
    plt.ylim(0, 1.2)
    plt.xlabel('$z$', fontsize = 13)
    plt.ylabel('Completeness', fontsize = 13)
    plt.plot(bin_x, compl[i], marker = '+', label = labels[i])

plt.title('AMICO-cosmoDC2')
plt.legend()
plt.savefig(outpath + 'z_' + 'nobckgcompl1D_plot_' +  matching + '.png', format='png', transparent=True)
        