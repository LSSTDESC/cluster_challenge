#!/usr/bin/env python
# coding: utf-8
# Study of DC2 stars to produce a mask for glaxies
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
from opening_catalogs_functions import *

print('Configuration arguments: ', str(sys.argv))
str_healpix_pixel = str(sys.argv[1])

###outpath
outpath = "/sps/lsst/users/tguillem/web/clusters/debug/DC2_mask/"
if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)

#get tract number
tract = sys.argv[1]

#function to print the content of hdf5 file
def get_all(name):
    print(name)

#FlexZBoost
#galaxy catalog
f1 = h5py.File(f'/sps/lsst/groups/desc/PZ/PhotoZDC2/run2.2i_dr6_v2/CATALOGS_W_TRUTH/Run2.2i_dr6_dereddened_tract_{str(tract)}_withtruez.hdf5', 'r')
f1.visit(get_all)
galaxy_id = f1['photometry/id']
ra = f1['photometry/ra']
dec = f1['photometry/dec']
redshift = f1['photometry/redshift']
mag_g = f1['photometry/mag_g_lsst']
mag_r = f1['photometry/mag_r_lsst']
mag_i = f1['photometry/mag_i_lsst']
mag_z = f1['photometry/mag_z_lsst']
mag_y = f1['photometry/mag_y_lsst']
match = f1["photometry/is_truth_match_good"]

#pdf catalog
f2 = h5py.File(f'/sps/lsst/groups/desc/PZ/PhotoZDC2/run2.2i_dr6_v2/FLEXZBOOST/PDF_FILES/FZBoost_tract_{str(tract)}_pdfs.hdf5', 'r')
f2.visit(get_all)
photoz_pdf = f2['data/yvals']
#print(photoz_pdf.shape)
#print(photoz_pdf[0])
xvals = f2['meta/xvals']
#print(xvals.shape)
#print(xvals[0])
zmode = f2['ancil/zmode']
#print(zmode.shape)
#print(zmode[0])

###recipe from D. Boutigny
dirDpdd =  '/sps/lsst/groups/desc/shared/DC2-prod/Run2.2i/dpdd/Run2.2i-dr6/object_dpdd_only/'
fileDpdd = dirDpdd +  f"object_dpdd_tract{str(tract)}.parquet"
columns = ["objectId", "clean", "good", "extendedness"]
f3 = pd.read_parquet(fileDpdd, columns=columns, engine='pyarrow')
clean = f3["clean"]
good = f3["good"]
extend = f3["extendedness"]
objectId = f3["objectId"]

#check that we have the exact same objects and in the same order in Flexzboost and in DR6 catalog
if (np.array(galaxy_id) - np.array(objectId)).sum() != 0 :
     print("Flexzboost and DR6 catalogs are different...")

galaxy_data_all = Table([galaxy_id, ra, dec, redshift, mag_g, mag_r, mag_i, mag_z, mag_y, zmode, photoz_pdf, match],
                        names=['id', 'ra', 'dec', 'redshift', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y', 'zmode', 'photoz_pdf', 'match'])

#define cut from the parquet tables
cut = clean & good & (np.array(extend).ravel()>0.5)
cut &= (np.array(mag_g).ravel()<27) & (np.array(mag_r).ravel()<27) & (np.array(mag_i).ravel()<27)
cut &= (np.array(mag_z).ravel()<27) & (np.array(mag_y).ravel()<27)
#cut &= np.array(redshift).ravel()<1.5

galaxy_data_cut = galaxy_data_all[np.array(cut).ravel()]
print(galaxy_data_cut[0:50])

table_xvals = Table([xvals,],names=('photoz_pdf_grid',))
#print(table_xvals)

galaxy_data_cut.write(outpath+'galaxies.fits')
table_xvals.write(outpath+'pdf_bins.fits')
sys.exit()


###
###galaxy_data_cut = galaxy_data_all[np.array(cut).ravel()]
###print(galaxy_data_cut[0:50])
###
###table_xvals = Table([xvals,],names=('photoz_pdf_grid',))
###print(table_xvals)
###
###galaxy_data_cut.write(outpath+'galaxies.fits')
###table_xvals.write(outpath+'pdf_bins.fits')
###sys.exit()             

#============HACK to study the mask
redshift=data['redshift']
mag_i=data['mag_i']
mag_i_photoz=data['mag_i_photoz']
mag_r_photoz=data['mag_r_photoz']
mag_z_photoz=data['mag_z_photoz']
mag_y_photoz=data['mag_y_photoz']
mag_z=data['mag_z']
mag_r=data['mag_r']
mag_g=data['mag_g']
mag_y=data['mag_y']
ra=data['ra']
dec=data['dec']
photoz_mask=data['photoz_mask']
galaxy_data_all = Table([ra, dec, redshift, mag_g, mag_i, mag_r, mag_z, mag_y, mag_i_photoz, mag_r_photoz, mag_z_photoz, mag_y_photoz, photoz_mask],names=('ra','dec', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y', 'mag_i_photoz', 'mag_r_photoz', 'mag_z_photoz', 'mag_y_photoz', 'photoz_mask'))
#galaxy_data_all = galaxy_data_all[galaxy_data_all['mag_i']>26.5]
#galaxy_data_all = galaxy_data_all[galaxy_data_all['photoz_mask']==False]
print(galaxy_data_all)
galaxy_data_all.write(outpath+'galaxies.fits')
sys.exit()
#============END HACK to study the mask

print('Loading catalog: DONE')
if photoz==True:
     #photoz grid
     photoz_pdf_bin_centers = Table(gc_truth.photoz_pdf_bin_centers)
     ###HACK because of the size issue in the photoz case
     redshift=data['redshift']
     photoz_mean=data['photoz_mean']
     photoz_mode=data['photoz_mode']
     photoz_mask=data['photoz_mask']
     print(photoz_mask.shape)
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
     mag_i_photoz=data['mag_i_photoz']
     mag_z_photoz=data['mag_z_photoz']
     mag_r_photoz=data['mag_r_photoz']
     mag_g_photoz=data['mag_g_photoz']
     mag_y_photoz=data['mag_y_photoz']
     ra=data['ra']
     dec=data['dec']
     #galaxy_id=data['galaxy_id']
     halo_mass=data['halo_mass']
     halo_id=data['halo_id']
     #apply mask
     #print(len(galaxy_id))
     #print(len(photoz_mask))
     #print(len(mag_i))
     #print(len(photoz_mean))
     #galaxy_id=galaxy_id[photoz_mask] ###already masked!
     ##ra=ra[photoz_mask]
     ##dec=dec[photoz_mask]
     ##halo_mass=halo_mass[photoz_mask]
     ##halo_id=halo_id[photoz_mask]
     ##redshift=redshift[photoz_mask]
     ##mag_g=mag_g[photoz_mask]
     ##mag_i=mag_i[photoz_mask]
     ##mag_r=mag_r[photoz_mask]
     ##mag_z=mag_z[photoz_mask]
     ##mag_y=mag_y[photoz_mask]
     ##mag_g_photoz=mag_g_photoz[photoz_mask]
     ##mag_i_photoz=mag_i_photoz[photoz_mask]
     ##mag_r_photoz=mag_r_photoz[photoz_mask]
     ##mag_z_photoz=mag_z_photoz[photoz_mask]
     ##mag_y_photoz=mag_y_photoz[photoz_mask]
     #galaxy_data_all = Table([galaxy_id, ra, dec, halo_mass, halo_id, redshift, mag_g, mag_i, mag_r, mag_z, mag_y, mag_g_photoz, mag_i_photoz, mag_r_photoz, mag_z_photoz, mag_y_photoz, photoz_mode, photoz_mean, photoz_pdf, photoz_median, photoz_odds, photoz_mode_ml_red_chi2, photoz_mode_ml],names=('galaxy_id','ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y',  'mag_g_photoz', 'mag_i_photoz', 'mag_r_photoz','mag_z_photoz', 'mag_y_photoz','photoz_mode', 'photoz_mean', 'photoz_pdf', 'photoz_median', 'photoz_odds', 'photoz_mode_ml_red_chi2', 'photoz_mode_ml'))
     galaxy_data_all = Table([ra, dec, halo_mass, halo_id, redshift, mag_g, mag_i, mag_r, mag_z, mag_y, mag_g_photoz, mag_i_photoz, mag_r_photoz, mag_z_photoz, mag_y_photoz, photoz_mode, photoz_mean, photoz_pdf, photoz_median, photoz_odds, photoz_mode_ml_red_chi2, photoz_mode_ml],names=('ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y',  'mag_g_photoz', 'mag_i_photoz', 'mag_r_photoz','mag_z_photoz', 'mag_y_photoz','photoz_mode', 'photoz_mean', 'photoz_pdf', 'photoz_median', 'photoz_odds', 'photoz_mode_ml_red_chi2', 'photoz_mode_ml'))
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

print(galaxy_data_all)
#conservative cuts to reduce size
#galaxy_data_all=galaxy_data_all[galaxy_data_all['redshift']<1.6]
#galaxy_data_all=galaxy_data_all[galaxy_data_all['mag_r']<27]

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
if DC2:
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
