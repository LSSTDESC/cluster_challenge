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

from opening_catalogs_functions import *

###clevar
import clevar
from clevar.catalog import ClCatalog
from clevar.match import ProximityMatch
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances

#outpath
##DC2_cat_name = 'cosmoDC2_v1.1.4'
DC2_cat_name = 'cosmoDC2_v1.1.4_small'

outpath = "/sps/lsst/users/tguillem/DESC/desc_may_2021/cluster_challenge/clevar_catalogs/cosmoDC2/" + DC2_cat_name + "/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

#richness and mass cuts
min_halo_mass = 10**14 #Msun
#gc_truth: catalog objects
truth_data, gc_truth = DC2_cat_open(DC2_cat_name, min_halo_mass, cluster_only=False)

#halo table
halo_data = truth_data[truth_data['is_central']==True]
galaxy_data = truth_data[truth_data['is_central']==False]

#c1 = ClCatalog('Cat1', ra=halo_data['ra'], dec=halo_data['dec'], z=halo_data['redshift'], mass=halo_data['halo_mass'], id=halo_data['halo_id'])
c1 = Table([halo_data['halo_id'],halo_data['ra'],halo_data['dec'],halo_data['redshift'],halo_data['halo_mass']],names=('id','ra','dec','z','mass'))

#c1.add_members(id=galaxy_data['galaxy_id'], id_cluster=galaxy_data['halo_id'], ra=galaxy_data['ra'], dec=galaxy_data['dec'])
c1_members = Table([galaxy_data['galaxy_id'],galaxy_data['halo_id'],galaxy_data['ra'],galaxy_data['dec'],galaxy_data['redshift']],names=('id','id_cluster','ra','dec','z')) 
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
