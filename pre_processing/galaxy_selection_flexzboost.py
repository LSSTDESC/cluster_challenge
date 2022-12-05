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
#print(ra.shape)
#print(ra[0])

#pdf catalog
f2 = h5py.File('/sps/lsst/users/boutigny/FLEXZBOOST/FZBoost_tract_3260_pdfs.hdf5', 'r')
f2.visit(get_all)
photoz_pdf = f2['data/yvals']
#print(yvals.shape)
#print(yvals[0])
xvals = f2['meta/xvals']
#print(xvals.shape)
#print(xvals[0])

galaxy_data_all = Table([galaxy_id, ra, dec, redshift, photoz_pdf], names=('galaxy_id','ra','dec', 'redshift', 'photoz_pdf'))
print(galaxy_data_all[0:4])

table_xvals = Table([xvals,],names=('photoz_pdf_grid',))
print(table_xvals)

galaxy_data_all.write('galaxies.fits')
table_xvals.write('pdf_bins.fits')

sys.exit()

#include or not photo-z pdfs
photoz = True

#gc and gc_truth: catalog objects / others are just tables
#halo_data, gc_truth = DC2_cat_open(DC2_cat_name, min_halo_mass, cluster_only=True)

#halo table
#mask = truth_data['is_central']==True
#mask = np.logical_and(truth_data['is_central']==True,truth_data['halo_id']==127200143167)
#halo_data = truth_data[mask]
#cosmoDC2_small: healpix_pixels = [9559,  9686,  9687,  9814,  9815,  9816,  9942,  9943, 10070, 10071, 10072, 10198, 10199, 10200, 10326, 10327, 10450]
#halo members
#galaxy_data_all = Table(gc_truth.get_quantities(['ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y', 'baseDC2/is_on_red_sequence_gr', 'baseDC2/is_on_red_sequence_ri', 'is_central', 'hostHaloMass',  'shear_1', 'shear_2',
#                                                 'convergence', 'baseDC2/sod_halo_mass'],
#                                                filters=['halo_mass > 10**13', 'mag_i<30', 'is_central==False']))# 'halo_id==127200143167']))

#galaxy_data_all = Table(gc_truth.get_quantities(['ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y']))
#galaxy_data_all = Table(gc_truth.get_quantities(['galaxy_id','ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y'],
#                                                filters=['is_central==False']))
#data = photoz_cat.get_quantities(['redshift', 'photoz_mode', 'photoz_mean','photoz_mask', 'mag_i', 'mag_z'], native_filters=partition_filter('healpix_pixel',['9559',  '9686',  '9687', '9814']))

print('Loading catalog')
#data = gc_truth.get_quantities(['galaxy_id','ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y','photoz_mode', 'photoz_mean', 'photoz_mask', 'photoz_pdf', 'photoz_median', 'photoz_odds', 'photoz_mode_ml_red_chi2', 'photoz_mode_ml'])#,native_filters=['healpix_pixel == 9559'])
#per healpix_pixel
print(str_healpix_pixel)
healpix_list = [str_healpix_pixel] 
data = gc_truth.get_quantities(['galaxy_id','ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y','photoz_mode', 'photoz_mean', 'photoz_mask', 'photoz_pdf', 'photoz_median', 'photoz_odds', 'photoz_mode_ml_red_chi2', 'photoz_mode_ml'], native_filters=['healpix_pixel == int(str_healpix_pixel)'])
###cosmoDC2 with photo-z
if photoz==True:
     data = gc_truth.get_quantities(['galaxy_id','ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y','photoz_mode', 'photoz_mean', 'photoz_mask', 'photoz_pdf', 'photoz_median', 'photoz_odds', 'photoz_mode_ml_red_chi2', 'photoz_mode_ml'], native_filters=partition_filter('healpix_pixel',healpix_list))
else:
     data = gc_truth.get_quantities(['galaxy_id','ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y'], native_filters=partition_filter('healpix_pixel',healpix_list))

###DC2 case
#cat = GCRCatalogs.load_catalog(DC2_cat_name)
#print(sorted(q for q in cat.list_all_quantities() if 'id' in q))
#data = cat.get_quantities(['objectId', 'ra', 'dec', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y','photoz_mode', 'photoz_mean', 'photoz_pdf', 'photoz_median', 'photoz_odds', 'photoz_mode_ml_red_chi2', 'photoz_mode_ml'], native_filters=['tract==4850'])
#data = cat.get_quantities(['objectId', 'ra', 'dec', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y','photoz_mode', 'photoz_mean', 'photoz_pdf', 'photoz_median', 'photoz_odds', 'photoz_mode_ml_red_chi2', 'photoz_mode_ml'], native_filters=tract_filter(healpix_list))
print('Loading catalog: DONE')
if photoz==True:
     #photoz grid
     photoz_pdf_bin_centers = Table(gc_truth.photoz_pdf_bin_centers)
     ###HACK because of the size issue in the photoz case
     redshift=data['redshift']
     photoz_mean=data['photoz_mean']
     photoz_mode=data['photoz_mode']
     photoz_mask=data['photoz_mask']
     #photoz_pdf=data['photoz_pdf']
     photoz_pdf=list(data['photoz_pdf'])
     photoz_median=data['photoz_median']
     photoz_odds=data['photoz_odds']
     photoz_mode_ml_red_chi2=data['photoz_mode_ml_red_chi2']
     photoz_mode_ml=data['photoz_mode_ml']
     mag_i=data['mag_i']
     mag_z=data['mag_z']
     mag_r=data['mag_r']
     mag_g=data['mag_g']
     mag_y=data['mag_y']
     ra=data['ra']
     dec=data['dec']
     galaxy_id=data['galaxy_id']
     halo_mass=data['halo_mass']
     halo_id=data['halo_id']
     #apply mask
     #print(len(galaxy_id))
     #print(len(photoz_mask))
     #print(len(mag_i))
     #print(len(photoz_mean))
     #galaxy_id=galaxy_id[photoz_mask] ###already masked!
     ra=ra[photoz_mask]
     dec=dec[photoz_mask]
     halo_mass=halo_mass[photoz_mask]
     halo_id=halo_id[photoz_mask]
     redshift=redshift[photoz_mask]
     mag_g=mag_g[photoz_mask]
     mag_i=mag_i[photoz_mask]
     mag_r=mag_r[photoz_mask]
     mag_z=mag_z[photoz_mask]
     mag_y=mag_y[photoz_mask]
     galaxy_data_all = Table([galaxy_id, ra, dec, halo_mass, halo_id, redshift, mag_g, mag_i, mag_r, mag_z, mag_y, photoz_mode, photoz_mean, photoz_pdf, photoz_median, photoz_odds, photoz_mode_ml_red_chi2, photoz_mode_ml],names=('galaxy_id','ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y','photoz_mode', 'photoz_mean', 'photoz_pdf', 'photoz_median', 'photoz_odds', 'photoz_mode_ml_red_chi2', 'photoz_mode_ml'))
else:
     photoz_pdf_bin_centers = Table(cat.photoz_pdf_bin_centers)
     photoz_mean=data['photoz_mean']
     photoz_mode=data['photoz_mode']
     photoz_pdf=list(data['photoz_pdf'])
     photoz_median=data['photoz_median']
     photoz_odds=data['photoz_odds']
     photoz_mode_ml_red_chi2=data['photoz_mode_ml_red_chi2']
     photoz_mode_ml=data['photoz_mode_ml']
     mag_i=data['mag_i']
     mag_z=data['mag_z']
     mag_r=data['mag_r']
     mag_g=data['mag_g']
     mag_y=data['mag_y']
     ra=data['ra']
     dec=data['dec']
     galaxy_id=data['objectId']
     ra=data['ra']
     dec=data['dec']
     mag_g=data['mag_g']
     mag_i=data['mag_i']
     mag_r=data['mag_r']
     mag_z=data['mag_z']
     mag_y=data['mag_y']
     galaxy_data_all = Table([galaxy_id, ra, dec, mag_g, mag_i, mag_r, mag_z, mag_y, photoz_mode, photoz_mean, photoz_pdf, photoz_median, photoz_odds, photoz_mode_ml_red_chi2, photoz_mode_ml],names=('galaxy_id','ra','dec', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y','photoz_mode', 'photoz_mean', 'photoz_pdf', 'photoz_median', 'photoz_odds', 'photoz_mode_ml_red_chi2', 'photoz_mode_ml'))

#conservative cuts to reduce size
#galaxy_data_all=galaxy_data_all[galaxy_data_all['redshift']<1.6]
galaxy_data_all=galaxy_data_all[galaxy_data_all['mag_r']<27]

###for DC2
#restrict to cosmoDC2 small footprint (common part)
#healpix selection
DC2 = False
if DC2:
     healpix_pixels = [9559,  9686,  9687,  9814,  9815,  9816,  9942,  9943, 10070, 10071, 10072, 10198, 10199, 10200, 10326, 10327, 10450]
     nside = 32
     filter_arr = []
     for element in galaxy_data_all:
          pix = hp.ang2pix(nside, element['ra'], element['dec'], lonlat=True)
          if ( pix in healpix_pixels ):
               filter_arr.append(True)
          else:
               filter_arr.append(False)
     galaxy_data_all = galaxy_data_all[filter_arr]
               
###redshift smearing
if not DC2:
     print('********Start redshift smearing********')
     redshift_smeared = []
     i=0
     for element in galaxy_data_all:
          i+=1
          sigma = 0.03*(1+element['redshift'])
          np.random.seed(i)
          smearing = np.random.normal(0, sigma, 1)
          z_smear = element['redshift'] + smearing
          redshift_smeared.append(z_smear)
          galaxy_data_all.add_column(redshift_smeared, index=3, name = 'redshift_smeared')
          print('********End redshift smearing********')

          ###mag selection
          print('********Start mag selection********')
          data_brightness_i = ascii.read("istar.asc")
          data_brightness_z = ascii.read("zstar.asc")
          filter_arr = []
          filter_arr2 = []
          for element in galaxy_data_all:
               index_i = mstar_i(element['redshift'])
               index_z = mstar_z(element['redshift'])
               if ( (element['mag_i'] < data_brightness_i[index_i][1]+2) and (element['mag_i']<25.47) ):
                    #if ( element['mag_i']<25.47 ):
                    filter_arr.append(True)
               else:
                    filter_arr.append(False)
          
               if ( (element['mag_z'] < data_brightness_z[index_z][1]+2) and (element['mag_z']<24.46) ):
               #if ( element['mag_z']<24.46 ):
                    filter_arr2.append(True)
               else:
                    filter_arr2.append(False)
          galaxy_data_all.add_column(filter_arr, name = 'mstar_i')          
          galaxy_data_all.add_column(filter_arr2, name = 'mstar_z')
          #galaxy_data_1 = galaxy_data_all[filter_arr]
          #galaxy_data_2 = galaxy_data_all[filter_arr2] 
          #print('********End mag selection********')

#print('********Start richness computation********')
#richness_i = []
#richness_z = []
#for halo in halo_data:
#     #_i
#     members_i = galaxy_data_1[galaxy_data_1['halo_id']==halo['halo_id']]
#     #print(members_i)
#     richness_i.append(len(members_i))
#     #_z
#     members_z = galaxy_data_2[galaxy_data_2['halo_id']==halo['halo_id']]
#     #print(members_z)
#     richness_z.append(len(members_z))
#     
#halo_data.add_column(richness_i, index=5, name = 'NGALS_i')
#halo_data.add_column(richness_z, index=5, name = 'NGALS_z')
#print(halo_data)
#print('********End richness computation********') 

#Save Tables
#halo_data.write(outpath+'halos.fits')
print(galaxy_data_all)
print(galaxy_data_all.info)
galaxy_data_all.write(outpath+'galaxies.fits')
photoz_pdf_bin_centers.write(outpath+'pdf_bins.fits')
#galaxy_data_1.write(outpath+'galaxy_data_1.fits')
#galaxy_data_2.write(outpath+'galaxy_data_2.fits')
#Check content
#print(halo_data)
#print(galaxy_data_all)
#print(photoz_pdf_bin_centers)
#mass richness (if sod_halo_mass available)
#plt.figure()
#plt.scatter(halo_data['NGALS_i'], halo_data['baseDC2/sod_halo_mass'], marker='.',color = 'blue', s=10, alpha=0.3, label='clusters')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim([10, 300])
#plt.ylim([2.0e13, 2.0e15])
#plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
#plt.savefig(outpath+'mass_richness.png', bbox_inches='tight')

sys.exit()
