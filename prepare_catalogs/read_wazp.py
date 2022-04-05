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
from clevar.catalog import ClCatalog
from clevar.match import ProximityMatch
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances

#outpath
algo = 'WaZP'
#6563 / 6684 / 6685 / 6688 / 6561
catalog_wazp = '6685'

outpath = "/sps/lsst/users/tguillem/DESC/desc_may_2021/cluster_challenge/clevar_catalogs/wazp/" + catalog_wazp + "/"

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
inpath = '/sps/lsst/users/tguillem/DESC/desc_may_2021/desc-data-portal/notebooks/dc2/'
wazp_cat_name = inpath + 'catalogs/' + catalog_wazp + '/wazp_cluster.fits'
wazp_members_cat_name = inpath + 'catalogs/' + catalog_wazp + '/wazp_membership.fits'
#truth matching
#matching = np.load('catalogs/wazp_truth_matching.npy')
#print(matching)
#test DES
#wazp_cat_name = 'catalogs/y1a1_dnf_wazp_v5.0.11.5239+47_clusters.fits'
#wazp_members_cat_name = 'catalogs/y1a1_dnf_wazp_v5.0.11.5239+47_members.fits'
wazp_data, wazp_members_data = wazp_cat_open(wazp_cat_name, wazp_members_cat_name, min_richness)
####General catalog properties
print("Number of WaZP clusters = ", len(wazp_data))
print("Number of WaZP cluster members = ", len(wazp_members_data))

#create clevar catalog
c1 = ClCatalog('Cat1', ra=wazp_data['ra'], dec=wazp_data['dec'], z=wazp_data['redshift'], mass = wazp_data['NGALS'], id=wazp_data['ID'])

#add members
c1.add_members(id=wazp_members_data['ID_g'], id_cluster=wazp_members_data['ID_CLUSTER'], ra=wazp_members_data['ra'], dec=wazp_members_data['dec'], pmem=wazp_members_data['PMEM'])

print(c1)
print(c1.members)
c1.write(outpath + 'ClCatalog.fits', overwrite=True)

sys.exit()
