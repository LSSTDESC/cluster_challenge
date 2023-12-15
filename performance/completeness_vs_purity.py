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

inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/amico_cosmoDC2/mag_i/w_rich/m14/cuts/'
sn_l = [3, 35, 4, 45, 5, 6, 7, 8]



# #---Select catalogs to match---##

outpath = '/pbs/home/n/namourou/test_jupyter/cluster_challenge/plots/'

print('Saving plots at :', outpath)
compl = np.empty(len(sn_l))
purity = np.empty(len(sn_l))
##########select case

nbins_x = 5
bin1 = np.linspace(0.0, 1.8, nbins_x) #For AMICO
for i in range(len(sn_l)):
    z_min = 0
    z_max = 1
    c1 = ClCatalog.read(inpath + 'c1_p_sn' + str(sn_l[i]) + '.fits', 'c1', full = True)
    c2 = ClCatalog.read(inpath + 'c2_p_sn' + str(sn_l[i]) + '.fits', 'c2', full = True)
    c_merged_12 = ClCatalog.read(inpath+'output_catalog_p_sn' + str(sn_l[i]) + '.fits', 'merged',  full=True)
    filter1 = np.logical_and(c_merged_12.data['cat2_z'] < z_max, c_merged_12.data['cat2_z'] > z_min)
    filter2 = np.logical_and(c2.data['z'] < z_max, c2.data['z'] > z_min)
    filter3 = np.logical_and(c1.data['z'] < z_max, c1.data['z'] > z_min)
    filter4 = np.logical_and(c_merged_12.data['cat1_z'] < z_max, c_merged_12.data['cat1_z'] > z_min)
    c_clusters = c1.data[filter3]
    c_halos_matched = c_merged_12[filter1]
    c_clusters_matched = c_merged_12[filter4]
    c_halos = c2.data[filter2]
    compl[i] = len(c_halos_matched)/len(c_halos)
    purity[i] = len(c_clusters_matched)/len(c_clusters)
plt.xlabel('completeness', fontsize = 13)
plt.ylabel('purity', fontsize = 13)
plt.xlim([0,1])
plt.ylim([0,1])
plt.plot(compl, purity, marker = '+')
plt.title('AMICO-cosmoDC2', fontsize = 13)
plt.legend()
plt.savefig(outpath + 'purity_vs_completeness_' + str(z_min) + '< z <' + str(z_max) + '.png')
