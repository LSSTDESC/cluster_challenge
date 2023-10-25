#!/usr/bin/env python
# coding: utf-8

###import
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
import numpy as np
from scipy.optimize import curve_fit
import sys
import os
import shutil
import pickle
import math
from scipy.optimize import minimize

###clevar
import clevar
from clevar.catalog import ClCatalog
from clevar.match import ProximityMatch
from clevar.match import get_matched_pairs
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances
from clevar.match_metrics.recovery import ClCatalogFuncs as r_cf
from clevar.match import output_matched_catalog

matching_folder = '/sps/lsst/groups/clusters/cluster_comparison_project/after_matching/cosmoDC2_wazp.cosmoDC2.fzb/v0_v0/member_matching/fshare_0.0_pref_shared_member_fraction/'
#define catalog names
catalog1 = 'wazp.cosmoDC2.fzb.v0.fits'
catalog2 = 'cosmoDC2.v0.fits'

outpath_base = '/pbs/home/t/tguillem/web/clusters/cluster_challenge/selection_function/completeness_versus_purity/'

#select the catalogs to match
wazp_cosmoDC2 = True
redmapper_cosmoDC2 = False
amico_cosmoDC2 = False
matching_selected = 'cross'

if wazp_cosmoDC2 == True:
     #matching_folder = matching_folder + 'wazp_cosmoDC2/'
     matching = 'WaZP cosmoDC2: NGALS>0, $m_{200c}>10^{13}$'
     outpath = outpath_base + 'wazp_cosmoDC2/'
elif redmapper_cosmoDC2 == True:
     #matching_folder = matching_folder + 'redmapper_cosmoDC2/'
     matching = 'redMaPPer cosmoDC2: $\Lambda>0$, $m_{200c}>10^{13}$'
     outpath = outpath_base + 'redmapper_cosmoDC2/'
elif amico_cosmoDC2 == True:
     matching_folder = matching_folder + 'amico_cosmoDC2/'
     matching = 'AMICO cosmoDC2: $m_{halo}>10^{13}$'
     outpath = outpath_base + 'amico_cosmoDC2/'
else:
     print('Catalog selection is wrong.')
     sys.exit()

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)

#load c1 and c2
c1 = ClCatalog.read_full(matching_folder + catalog1)
c2 = ClCatalog.read_full(matching_folder + catalog2)

#create a merged catalog for the cross-matched pairs
output_matched_catalog(matching_folder+catalog1, matching_folder+catalog2,matching_folder+'output_catalog_12.fits', c1, c2, matching_type='cross', overwrite=True)
c_merged_12 = ClCatalog.read(matching_folder+'output_catalog_12.fits', name='NameOfCatalog', full=True)

#update plot title
matching = 'WaZP cosmoDC2: $m_{200c}>10^{14}$, z<1.15, NGALS>0'
labels=['z incl.']
colors=['black','red','blue','purple']

#completeness versus SNR
#nbins_x = 8
#snrbins = [3,4,5,7,9,12,15,20,30]
nbins_x = 7
snrbins = [3,4,5,7,9,12,15,20]
compl_snr = np.empty(nbins_x)

#prepare tables
c_merged=c_merged_12[c_merged_12.data['cat2_log10_mass']<100]
c_merged=c_merged[c_merged['cat2_z_cl']<1.15]
c_merged=c_merged[c_merged['cat1_z_cl']<1.15]
c_merged_cut=c_merged[c_merged['cat2_log10_mass']>14]
                      
c_halos = c2.data[c2.data['log10_mass']>13]
c_halos = c_halos[c_halos['z_cl']<1.15]
c_halos_cut = c_halos[c_halos['log10_mass']>14]

c_clusters = c1.data[c1.data['snr_cl']>0]
c_clusters = c_clusters[c_clusters['z_cl']<1.15]

for i in range(0,nbins_x):
     print('-----'+str(i))
     cut1 = snrbins[i]
     cut2 = snrbins[i+1]
     c_halos_matched = c_merged_cut[c_merged_cut['cat1_snr_cl']>cut1]
     #c_halos = c2.data
     n_halos_matched = len(c_halos_matched)
     n_halos = len(c_halos_cut)
     print(n_halos_matched)
     print(n_halos)
     compl_snr[i] = round(n_halos_matched/n_halos,4)
     print(compl_snr[i])

bin_x_snr = np.empty(nbins_x)
for ix in range(nbins_x):
     bin_x_snr[ix] = 0.5 * (snrbins[ix] + snrbins[ix+1])
     
plt.figure()
plt.scatter(bin_x_snr, compl_snr, label=labels[0], color=colors[0], marker= ".", s=30)
plt.plot(bin_x_snr, compl_snr, color=colors[0])
plt.ylim(0, 1.2)
plt.xlim(2,20)
plt.xlabel('SNR')
plt.ylabel('Completeness')
plt.title(matching)
#plt.legend()
plt.savefig(outpath+"completeness_snr.png", bbox_inches='tight')
plt.close()

#purity versus SNR
#nbins_x = 8
#snrbins = [3,4,5,6,7,8,9,10,11]
purity_snr = np.empty(nbins_x)
for i in range(0,nbins_x):
     print('-----'+str(i))
     cut1 = snrbins[i]
     cut2 = snrbins[i+1]
     c_clusters_matched = c_merged[c_merged['cat1_snr_cl']>cut1]
     c_clusters_all = c_clusters[c_clusters['snr_cl']>cut1]
     n_clusters_matched = len(c_clusters_matched)
     n_clusters_all = len(c_clusters_all)
     print(n_clusters_matched)
     print(n_clusters_all)
     purity_snr[i] = round(n_clusters_matched/n_clusters_all,4)
     print(purity_snr[i])

bin_x_snr = np.empty(nbins_x)
for ix in range(nbins_x):
     bin_x_snr[ix] = 0.5 * (snrbins[ix] + snrbins[ix+1])
     
plt.figure()
plt.scatter(bin_x_snr, purity_snr, label=labels[0], color=colors[0], marker= ".", s=30)
plt.plot(bin_x_snr, purity_snr, color=colors[0])
plt.ylim(0, 1.2)
plt.xlim(2,20)
plt.xlabel('SNR')
plt.ylabel('Purity')
plt.title(matching)
#plt.legend()
plt.savefig(outpath+"purity_snr.png", bbox_inches='tight')
plt.close()

#completeness versus purity
plt.figure()
plt.scatter(compl_snr, purity_snr, label=labels[0], color=colors[0], marker= ".", s=30)
plt.plot(compl_snr, purity_snr, color=colors[0])
plt.ylim(0.3, 1)
plt.xlim(0, 1)
plt.xlabel('Completeness')
plt.ylabel('Purity')
plt.title(matching)
#plt.legend()
plt.savefig(outpath+"purity_versus_completeness.png", bbox_inches='tight')
plt.close()

print('DONE')
sys.exit()
