#!/usr/bin/env python
# coding: utf-8

###import
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
#import esutil
import sys
import os
import shutil
import pickle
import h5py
import pandas as pd

###clevar
import clevar
from clevar.catalog import ClCatalog
from clevar.match import ProximityMatch
from clevar.match import get_matched_pairs
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances
from clevar.match_metrics.scaling import ClCatalogFuncs as s_cf
from clevar.match import output_matched_catalog

###plot style
plt.rcParams['figure.figsize'] = [9.5, 6]
plt.rcParams.update({'font.size': 18})
#plt.rcParams['figure.figsize'] = [10, 8] for big figures
#########################################

outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/selection_function/halos/"
if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)

#inpath = "/sps/lsst/users/ccombet/SkySim5000/hdf5/"
inpath = "/sps/lsst/users/tguillem/catalogs/SkySim5000/hdf5/"

#read hdf5 files
with pd.HDFStore(os.path.join(inpath,f'skysim_halos_z=0-1.20_mfof_gt_1.00e+13_small.hdf5')) as store:
     halo_data = store['skysim']
     halo_metadata = store.get_storer('skysim').attrs.metadata
#print(halo_data['baseDC2/sod_halo_mass'])

#rename
halo_data.rename(columns={'baseDC2/sod_halo_mass': 'M200c', 'richness': 'NGALS', 'richness_i': 'NGALS_i', 'richness_z': 'NGALS_z'}, inplace=True)
print(halo_data)
#fix M200c
halo_data['M200c'] = halo_data['M200c']/0.71
#mass richness
thislist = ["NGALS", "NGALS_i", "NGALS_z"]
for richness in thislist:
     plt.figure()
     plt.scatter(halo_data[richness], halo_data['M200c'], marker='.',color = 'blue', s=10, alpha=0.3, label='clusters')
     plt.xscale('log')
     plt.yscale('log')
     plt.xlim([1, 100])
     plt.ylim([1.0e12, 2.0e15])
     plt.xlabel(richness)
     plt.ylabel('M200c')
     plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
     plt.savefig(outpath+'mass_'+richness+'.png', bbox_inches='tight') 
     plt.close()
     #richness versus redshift
     plt.figure()
     plt.scatter(halo_data['redshift'], halo_data[richness], marker='.',color = 'blue', s=10, alpha=0.3, label='clusters')
     #plt.xscale('log')
     plt.yscale('log')
     plt.xlim([0, 1.3])
     plt.ylim([1, 100])
     plt.xlabel('redshift')
     plt.ylabel(richness)
     plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
     plt.savefig(outpath+''+richness+'_redshift.png', bbox_inches='tight') 
     plt.close()

#richness plots in (m,z) bins
zbins = [0,0.5,0.75,1.0,1.2]
ybins = 10**np.linspace(13, 14.5, 9)
n_z = len(zbins)-1
n_y = len(ybins)-1
print(ybins)
for j in range(0,n_y):
     print(ybins[j])
     
for i in range(0,n_z):
     cut1 = zbins[i]
     cut2 = zbins[i+1]
     filter1 = np.logical_and(halo_data['redshift'] > cut1, halo_data['redshift'] < cut2)
     halos_1 = halo_data[filter1]
     for j in range(0,n_y):
          cut3 = ybins[j]
          cut4 = ybins[j+1]
          print(cut3)
          filter2 = np.logical_and(halos_1['M200c'] > cut3, halos_1['M200c'] < cut4)
          halos = halos_1[filter2]
          print(halos)
          #richness
          nbins = 30
          bin_range = [0,150]
          plt.figure()
          plt.hist(halos['NGALS'], range=bin_range, bins=nbins, label='NGALS', histtype='step', color = 'black', density=True, stacked=True)
          #plt.hist(halos['NGALS_i'], range=bin_range, bins=nbins, label='NGALS_i', histtype='step', color = 'red')
          #plt.hist(halos['NGALS_z'], range=bin_range, bins=nbins, label='NGALS_z', histtype='step', color = 'blue')
          plt.xlabel("NGALS");
          plt.ylabel("halos")
          #plt.xscale('log')
          #plt.yscale('log')
          plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
          plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
          #plt.title('z '+str(cut1)+'-'+str(cut2)+'mass '+str(cut1)+'-'+str(cut2))
          plt.legend(title = '', loc='upper right')
          plt.savefig(outpath+'richness_redshift_bin_'+str(i)+'_mass_bin_'+str(j)+'.png')
          plt.close()
                                                                                          
sys.exit()
     
#M200c vs halo_mas
plt.figure()
plt.scatter(halo_data['halo_mass'], halo_data['M200c'], marker='.',color = 'blue', s=10, alpha=0.3, label='clusters')
plt.xscale('log')
plt.yscale('log')
plt.xlim([1.0e13, 1.0e15])
plt.ylim([1.0e12, 2.0e15])
plt.xlabel('halo_mass')
plt.ylabel('M200c')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'mass_mass.png', bbox_inches='tight') 

zbins = [0,0.5,0.75,1.0,1.25,1.5]
for i in range(0,5):
     cut1 = zbins[i]
     cut2 = zbins[i+1]
     filter1 = np.logical_and(halo_data['redshift'] > cut1, halo_data['redshift'] < cut2)
     halos = halo_data[filter1]
     #filter2 = np.logical_and(filter1,halo_data['NGALS_i']>25)
     #halos = halo_data[filter2]
     
     #2d histo from numpy
     ybins = 10**np.linspace(13, 15, 20)
     xbins = [25,30,35,40,50,60,100]
     #M200c
     counts, xedges, yedges = np.histogram2d(halos['NGALS_i'], halos['M200c'], bins=(xbins, ybins))
     fig, ax = plt.subplots()
     ax.pcolormesh(xbins, ybins, counts.T, cmap='jet')
     ax.set_yscale('log')
     ax.set_xscale('log')
     ax.set_xlabel('NGALS_i')
     ax.set_ylabel('M200c')
     fig.savefig(outpath+'mass_richness_2D_z_bin_'+str(i)+ '.png', bbox_inches='tight')

     file = open(outpath+'richness_z_bin_'+str(i)+ '.p', "wb")
     pickle.dump(counts, file)
     pickle.dump(xedges, file)
     pickle.dump(yedges, file)
     file.close() 

     #halo_mass
     counts, xedges, yedges = np.histogram2d(halos['NGALS_i'], halos['halo_mass'], bins=(xbins, ybins))
     fig, ax = plt.subplots()
     ax.pcolormesh(xbins, ybins, counts.T, cmap='jet')
     ax.set_yscale('log')
     ax.set_xscale('log')
     ax.set_xlabel('NGALS_i')
     ax.set_ylabel('m_halo')
     fig.savefig(outpath+'mass_richness_halo_mass_2D_z_bin_'+str(i)+ '.png', bbox_inches='tight')

     file = open(outpath+'richness_halo_mass_z_bin_'+str(i)+ '.p', "wb")
     pickle.dump(counts, file)
     pickle.dump(xedges, file)
     pickle.dump(yedges, file)
     file.close()

     #halo_mass vs M200c scatter
     plt.figure()
     plt.scatter(halos['halo_mass'], halos['M200c'], marker='.',color = 'blue', s=10, alpha=0.3, label='clusters')
     plt.xscale('log')
     plt.yscale('log')
     plt.xlim([5.0e13, 1.0e15])
     plt.ylim([1.0e13, 1.0e15])
     plt.xlabel('halo_mass')
     plt.ylabel('M200c')
     plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
     plt.title('z '+str(cut1)+'-'+str(cut2))
     plt.savefig(outpath+'mass_mass_z_bin_'+str(i)+'.png', bbox_inches='tight') 

     #ratio plots
     plt.figure()
     bin_range = [0.5,1.3]
     nbins = 80
     plt.hist(halos['M200c']/halos['halo_mass'], range=bin_range, bins=nbins, histtype='step', color = 'black')
     plt.xlabel("M200c/halo_mass");
     plt.ylabel("halos")
     plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
     plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
     #plt.legend(title = '', loc='upper right')
     plt.title('z '+str(cut1)+'-'+str(cut2)) 
     plt.savefig(outpath+'mass_ratio_z_bin_'+str(i)+'.png')
     plt.close()

     #richness
     nbins = 20
     bin_range = [0,100]
     plt.figure()
     plt.hist(halos['NGALS'], range=bin_range, bins=nbins, label='NGALS', histtype='step', color = 'black')
     plt.hist(halos['NGALS_i'], range=bin_range, bins=nbins, label='NGALS_i', histtype='step', color = 'red')
     plt.hist(halos['NGALS_z'], range=bin_range, bins=nbins, label='NGALS_z', histtype='step', color = 'blue')
     plt.xlabel("NGALS");
     plt.ylabel("halos")
     #plt.xscale('log')
     #plt.yscale('log')
     plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
     plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
     plt.title('z '+str(cut1)+'-'+str(cut2))
     plt.legend(title = '', loc='upper right')
     plt.savefig(outpath+'richness_bin_'+str(i)+".png")
     plt.close()



     
sys.exit()
