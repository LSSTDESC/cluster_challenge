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
import math

from opening_catalogs_functions import *

###clevar
import clevar
from clevar import ClCatalog
from clevar.match import ProximityMatch
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances

#outpath
DC2_cat_name = 'cosmoDC2_v1.1.4_image'
#DC2_cat_name = 'cosmoDC2_v1.1.4_small'

outpath = "/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/cosmoDC2/m200c_gt_13.0/" + DC2_cat_name + "/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

##########to add m200c
#method 1
halo_data = Table.read('/sps/lsst/groups/clusters/dc2/cosmoDC2_v1.1.4/extragal/halos/halos_m200c_13.0.fits')
members_m200c = Table.read('/sps/lsst/groups/clusters/dc2/cosmoDC2_v1.1.4/extragal/halos/halos_m200c_13.0_members.fits')
#halos_m200c = Table.read('/sps/lsst/users/maguena/cats/dc2/cosmoDC2_v1.1.4/extragal/halos_full/halos_m12.885.fits')
#halo_data = Table.read('/sps/lsst/users/maguena/cats/dc2/cosmoDC2_v1.1.4/extragal/halos_full/halos_m12.885_mtskysim.fits')
print(members_m200c.colnames)
#print(halos_m200c)
#sys.exit()

#method 2
#read hdf5 files
#inpath = "/sps/lsst/users/ccombet/SkySim5000/hdf5/"
#inpath = "/tmp/tguillem/"
#with pd.HDFStore(os.path.join(inpath,f'skysim_halos_z=0-1.20_mfof_gt_1.00e+13_small.hdf5')) as store:
#     halos_m200c = store['skysim']
#     halo_metadata = store.get_storer('skysim').attrs.metadata
#     #print(halo_data['baseDC2/sod_halo_mass'])
     
#rename
#halos_m200c.rename(columns={'baseDC2/sod_halo_mass': 'm200c', 'richness': 'NGALS', 'richness_i': 'NGALS_i', 'richness_z': 'NGALS_z'}, inplace=True)
#print(halos_m200c)
#fix M200c
#halos_m200c['m200c'] = halos_m200c['m200c']/0.71

#method 3
#read from GCR
#DC2_cat_name2 = 'skysim5000_v1.1.1'
#min_halo_mass = 10**13 #Msun
#gc_truth2 = GCRCatalogs.load_catalog(DC2_cat_name2)
#query = GCRCatalogs.GCRQuery('(is_central == True) & (halo_mass > ' + str(min_halo_mass) +')')
#skysim_data = Table(gc_truth2.get_quantities(['redshift','halo_mass','halo_id','galaxy_id','ra','dec','is_central','baseDC2/sod_halo_mass'],[query]))
#c0 = Table([skysim_data['halo_id'],skysim_data['ra'],skysim_data['dec'],skysim_data['redshift'],skysim_data['baseDC2/sod_halo_mass']],names=('id','ra','dec','z','m200c'))
#fix M200c
#c0['m200c'] = c0['m200c']/0.71
##########

#mass cut
min_halo_mass = 10**12.885 #Msun
#gc_truth: catalog object
#truth_data, gc_truth = DC2_cat_open(DC2_cat_name, min_halo_mass, cluster_only=False)
#galaxy_data = truth_data[truth_data['is_central']==False]

c1 = Table([halo_data['halo_id'],halo_data['ra'],halo_data['dec'],halo_data['redshift_true'],halo_data['mass_fof'],halo_data['m200c']],names=('id','ra','dec','z','mass','m200c'))
#add log(m) column
c1.add_column(1.0, name='log_mass', index=5)
for i in range(0,len(c1)):
     c1['log_mass'][i]=np.log10(c1['mass'][i])
#add log(m200c) column
c1.add_column(1.0, name='log_m200c', index=7)
for i in range(0,len(c1)):
     c1['log_m200c'][i]=np.log10(c1['m200c'][i])

#c1.add_members(id=galaxy_data['galaxy_id'], id_cluster=galaxy_data['halo_id'], ra=galaxy_data['ra'], dec=galaxy_data['dec'])
#from GCR
#c1_members = Table([galaxy_data['galaxy_id'],galaxy_data['halo_id'],galaxy_data['ra'],galaxy_data['dec'],galaxy_data['redshift']],names=('id','id_cluster','ra','dec','z')) 
#c1_members.add_column(1.0, name='pmem', index=5)
#from already-prepared table
c1_members = Table([members_m200c['galaxyID'],members_m200c['halo_id'],members_m200c['ra'],members_m200c['dec'],members_m200c['redshift_true'],members_m200c['mag_true_g'],members_m200c['mag_true_r'],members_m200c['mag_true_i'],members_m200c['mag_true_z'],members_m200c['mag_true_y']],names=('id','id_cluster','ra','dec','z','mag_g','mag_r','mag_i','mag_z','mag_y'))
#members_m200c['pmem_nfw2d']
c1_members.add_column(1.0, name='pmem', index=5)

#print(c1)
#print(c1.members)
#c1.write(outpath + 'ClCatalog.fits', overwrite=True)
#c1.members.write(outpath + 'ClCatalog_members.fits', overwrite=True) 
print(c1)
print(c1_members)
c1.write(outpath + 'Catalog.fits', overwrite=True)
c1_members.write(outpath + 'Catalog_members.fits', overwrite=True)

sys.exit()
