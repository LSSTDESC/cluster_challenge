#!/usr/bin/env python
# coding: utf-8

###import
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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

from wazp_functions import *

###clevar
import clevar
from clevar import ClCatalog
from clevar.match import ProximityMatch
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances

#outpath
algo = 'WaZP'
#6563 / 6684 / 6685 / 6688 / 6561 / 
catalog_wazp = '6980'

outpath = "/sps/lsst/groups/clusters/cluster_comparison_project/debug/wazp/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

###WaZP specific format
#cat_wazp = fits.open('catalogs/wazp_cluster.fits')
#cat_wazp_table = cat_wazp[0]
#cat_wazp_original = Table.read('catalogs/wazp_cluster.fits')
#cat_wazp = 
#print(cat_wazp.info)
#richness and mass cuts
min_richness = 0
#inpath = '/sps/lsst/users/tguillem/DESC/desc_may_2021/desc-data-portal/notebooks/dc2/'
#wazp_cat_name = inpath + 'catalogs/wazp_full/' + catalog_wazp + '/wazp_cluster.fits'
#wazp_members_cat_name = inpath + 'catalogs/wazp_full/' + catalog_wazp + '/wazp_membership.fits'
#NEW: 07/04/23
#inpath = '/sps/lsst/groups/clusters/wazp_validation_project/6948/'
#cosmoDC2 flexzboost
inpath = '/sps/lsst/users/tguillem/DESC/desc_may_2021/desc-data-portal/notebooks/dc2/catalogs/wazp_full/6980/'
wazp_cat_name = inpath + 'wazp_cluster.fits'
wazp_members_cat_name = inpath + 'wazp_membership.fits'
#truth matching
#matching = np.load('catalogs/wazp_truth_matching.npy')
#print(matching)
#test DES
#wazp_cat_name = 'catalogs/y1a1_dnf_wazp_v5.0.11.5239+47_clusters.fits'
#wazp_members_cat_name = 'catalogs/y1a1_dnf_wazp_v5.0.11.5239+47_members.fits'

###method 1: open local .fits files
wazp_data, wazp_members_data = wazp_cat_open(wazp_cat_name, wazp_members_cat_name, min_richness)

###method 2: open GCR catalogs
#print('\n'.join(sorted(c for c in GCRCatalogs.get_available_catalogs(include_default_only=False) if c.startswith('cosmoDC2'))))
#wazp_cat_name = 'cosmoDC2_v1.1.4_wazp_v1.0_flexzboost_v1'
#wazp_data, wazp_members_data, gc = wazp_cat_open_gcr(wazp_cat_name, min_richness) 

#wazp_data=wazp_data[wazp_data['NGALS']>5]
####General catalog properties
print("Number of WaZP clusters = ", len(wazp_data))
print("Number of WaZP cluster members = ", len(wazp_members_data))
#wazp_data=wazp_data[0:3]
#wazp_data.pprint_all()
#sys.exit()

#create clevar catalog
#c1 = ClCatalog('Cat1', ra=wazp_data['ra'], dec=wazp_data['dec'], z=wazp_data['redshift'], mass = wazp_data['NGALS'], id=wazp_data['ID'])
c1 = Table([wazp_data['ID'],wazp_data['ra'],wazp_data['dec'],wazp_data['redshift'],wazp_data['NGALS'],wazp_data['SNR']],names=('id','ra','dec','z','mass','snr'))

#add log(mass) column
c1.add_column(1.0, name='log_mass', index=6)
for i in range(0,len(c1)):
          c1['log_mass'][i]=np.log10(c1['mass'][i])

#add members
#c1.add_members(id=wazp_members_data['ID_g'], id_cluster=wazp_members_data['ID_CLUSTER'], ra=wazp_members_data['ra'], dec=wazp_members_data['dec'], pmem=wazp_members_data['PMEM'])
c1_members = Table([wazp_members_data['ID_g'],wazp_members_data['ID_CLUSTER'],wazp_members_data['ra'],wazp_members_data['dec'],wazp_members_data['redshift'],wazp_members_data['PMEM'],wazp_members_data['mag_r'],wazp_members_data['mag_i'],wazp_members_data['mag_z']],names=('id','id_cluster','ra','dec','z','pmem','mag_r','mag_i','mag_z'))

#add m* columns
data_brightness_i = ascii.read("mstar_files/istar.asc")
data_brightness_z = ascii.read("mstar_files/zstar.asc")
filter_arr = []
filter_arr2 = []
for element in c1_members:
     index_i = mstar_i(element['z'])
     index_z = mstar_z(element['z'])
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
c1_members.add_column(filter_arr, name = 'mstar_i')
c1_members.add_column(filter_arr2, name = 'mstar_z')

#print(c1)
#print(c1.members)
#c1.write(outpath + 'ClCatalog.fits', overwrite=True)
#c1.members.write(outpath + 'ClCatalog_members.fits', overwrite=True)
print(c1)
print(c1_members)
c1.write(outpath + 'Catalog.fits', overwrite=True)
c1_members.write(outpath + 'Catalog_members.fits', overwrite=True)

sys.exit()
