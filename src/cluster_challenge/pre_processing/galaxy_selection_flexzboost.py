#!/usr/bin/env python
# coding: utf-8
# Selection of galaxies from GCR catalogs to produce .FITS files per healpix pixel as inputs to galaxy cluster algorithms
# Author: Thibault Guillemin

###import
import GCRCatalogs
from GCRCatalogs.helpers.base_filters import sample_filter, partition_filter
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
#from GCRCatalogs.helpers.base_filters import partition_filter
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
import sys
import os
import shutil
import pickle
import healpy as hp
import h5py
import pandas as pd 

sys.path.append('../prepare_catalogs/')
from wazp_functions import *

###cluster_validation
#from cluster_validation.opening_catalogs_functions import *
#from cluster_validation.wazp_functions import *

print('Configuration arguments: ', str(sys.argv))
str_healpix_pixel = str(sys.argv[1])

###outpath
#outpath = "/sps/lsst/users/tguillem/DESC/desc_may_2021/desc-data-portal/notebooks/dc2/tables/cosmoDC2_small_photozs_flexzboost/"+str_healpix_pixel+"/" 
outpath = "/sps/lsst/users/tguillem/DESC/desc_may_2021/desc-data-portal/notebooks/dc2/tables/DC2/"+str_healpix_pixel+"/"
if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)

#mu, sigma = 0, 1 # mean and standard deviation
#np.random.seed(2022)
#smearing = np.random.normal(mu, sigma, 10000000)#10.000.000

###catalogs
#DC2_cat_name = 'cosmoDC2_v1.1.4'
#DC2_cat_name = 'cosmoDC2_v1.1.4_small'
#DC2_cat_name = 'cosmoDC2_v1.1.4_small_with_photozs_v1'
#DC2_cat_name = 'cosmoDC2_v1.1.4_small_with_photozs_flexzboost_v1'
#DC2_cat_name = 'dc2_object_run2.2i_dr6_v2_with_addons'
#DC2_cat_name = 'skysim5000_v1.1.1_small'
#DC2_cat_name = 'skysim5000_v1.1.1'
#richness and mass cuts
min_halo_mass = 10**13#Msun

def get_all(name):
        print(name)

#FlexZBoost
#galaxy catalog
f1 = h5py.File('/sps/lsst/users/boutigny/FLEXZBOOST/Run2.2i_dr6_dereddened_tract_3260_withtruez.hdf5', 'r')
f1.visit(get_all)
galaxy_id = f1['photometry/id']
ra = f1['photometry/ra']
dec = f1['photometry/dec']
redshift = f1['photometry/redshift']
mag_i = f1['photometry/mag_i_lsst']
mag_z = f1['photometry/mag_z_lsst']
print(ra.shape)
print(ra[0])

#pdf catalog
f2 = h5py.File('/sps/lsst/users/boutigny/FLEXZBOOST/FZBoost_tract_3260_pdfs.hdf5', 'r')
f2.visit(get_all)
photoz_pdf = f2['data/yvals']
print(photoz_pdf.shape)
print(photoz_pdf[0])
xvals = f2['meta/xvals']
print(xvals.shape)
print(xvals[0])
zmode = f2['ancil/zmode']
print(zmode.shape)
print(zmode[0])

galaxy_data_all = Table([galaxy_id, ra, dec, redshift, photoz_pdf, mag_i, mag_z], names=('galaxy_id','ra','dec', 'redshift', 'photoz_pdf', 'mag_i', 'mag_z'))
print(galaxy_data_all[0:50])

table_xvals = Table([xvals,],names=('photoz_pdf_grid',))
print(table_xvals)

#galaxy_data_all.write('galaxies.fits')
#table_xvals.write('pdf_bins.fits')

#########validation plots
###Photoz plot
#galaxy_data_all = galaxy_data_all[galaxy_data_all['mag_i']<25.5]
outpath = '/pbs/home/t/tguillem/web/clusters/cluster_challenge/debug/'
#plt.figure()
#plt.hist(galaxy_data_all['redshift'], density=False,range = (0, 2), bins = 40, label="redshift (mag_i<25.5)", histtype='step', color = 'black');
#plt.hist(galaxy_data_1['photoz_mean'], density=False, range = (0, 2), bins = 40, label="pz_mean (mag_i<25.5)", histtype='step', color = 'red');
#plt.hist(galaxy_data_1['photoz_mode'], density=False, range = (0, 2), bins = 40, label="pz_mode (mag_i<25.5)", histtype='step', color = 'blue');
#plt.hist(galaxy_data_2['redshift'], density=False, range = (0, 2), bins = 40, label="pz (mag_z<24.5)", histtype='step', color = 'blue');
#plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
# Customize the minor grid
#plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.5', color='grey')
#plt.xlabel("z");
#plt.ylabel("galaxies / 0.05 dz");
#plt.legend([g1,g2], ["dr3", "cosmoDC2"], title = 'Photo-z', loc='upper right')
#plt.legend(title = 'DR6 Flexzboost (tract 3260)', loc='upper right')
#outpath = "/sps/lsst/users/tguillem/DESC/desc_may_2021/desc-data-portal/notebooks/dc2/plots/photoz/bpz/"
#sys.mkdir(outpath);
#plt.savefig(outpath+"photoz.png", bbox_inches='tight')
#plt.close()

#good galaxy selection
#print('********Start mag selection********')
#data_brightness_i = ascii.read("istar.asc")
data_brightness_z = ascii.read("/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/prepare_catalogs/mstar_files/zstar.asc")
filter_arr = []
for element in galaxy_data_all:
     index_z = mstar_z(element['redshift'])
     if ( (element['mag_z'] < data_brightness_z[index_z][1]+2) and (element['mag_z']<24.46)):
          filter_arr.append(True)
     else:
          filter_arr.append(False)
galaxy_data_all = galaxy_data_all[filter_arr]
#zmode = zmode.flatten()
#zmode = zmode[filter_arr]
zmode = zmode[:,0]
zmode = zmode[filter_arr]

#pdf shape
#select one galaxy
#galaxy_sel = galaxy_data_all[galaxy_data_all['galaxy_id']==14337631626200577]
for i_gal in range(100):
     galaxy_sel = galaxy_data_all[i_gal]
     #print(galaxy_sel)

     i_gal_str = str(i_gal)

     #prepare pdf for plot
     pdf = 0.01*galaxy_sel['photoz_pdf'].flatten()

     #plot
     plt.figure()
     plt.hist(galaxy_sel['redshift'], density=False, range = (0, 3), bins = 300, label="redshift", histtype='step', color = 'red')
     plt.hist(zmode[i_gal], density=False, range = (0, 3), bins = 300, label="zmode", histtype='step', color = 'blue')
     #plt.hist(pdf, density=False,range = (0, 3), bins = 300, label="photoz_pdf", histtype='step', color = 'black')
     #plt.hist(bins_x, bins_x, weights=pdf, label="photoz_pdf", histtype='step', color = 'black')
     plt.scatter(xvals[0], pdf, label="photoz_pdf", color= "black", marker= ".", s=30)
     #table_xvals
     plt.ylim(0,0.2)
     x_av = zmode[i_gal]
     plt.xlim(x_av-0.5,x_av+0.5)
     plt.legend()
     plt.savefig(outpath+"pdf_gal_"+i_gal_str+".png", bbox_inches='tight')
     plt.close()

#correlation plots
plt.figure()
plt.scatter(galaxy_data_all['redshift'],zmode, marker='.',color = 'blue', s=0.5, alpha=0.3, label='galaxies')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xlabel('true z')
plt.ylabel('zmode')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'redshift_correlation.png', bbox_inches='tight')

plt.figure()
plt.scatter(galaxy_data_all['redshift'],zmode-galaxy_data_all['redshift'], marker='.',color = 'blue', s=0.5, alpha=0.3, label='galaxies')
plt.xlim([0, 1.5])
plt.ylim([-0.3, 0.3])
plt.xlabel('zmode')
plt.ylabel('zmode-true_z')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'redshift_difference.png', bbox_inches='tight')

plt.figure()
plt.scatter(zmode,zmode-galaxy_data_all['redshift'], marker='.',color = 'blue', s=1, alpha=0.3, label='galaxies')
plt.xlim([0.1, 0.5])
plt.ylim([-0.1, 0.1])
plt.xlabel('zmode')
plt.ylabel('zmode-true_z')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'redshift_difference_low.png', bbox_inches='tight')

plt.figure()
plt.scatter(zmode,zmode-galaxy_data_all['redshift'], marker='.',color = 'blue', s=0.5, alpha=0.3, label='galaxies')
plt.xlim([0.8, 1.5])
plt.ylim([-0.3, 0.3])
plt.xlabel('zmode')
plt.ylabel('zmode-true_z')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'redshift_difference_high.png', bbox_inches='tight')

sys.exit()
