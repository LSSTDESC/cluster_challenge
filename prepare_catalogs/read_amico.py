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
algo = 'AMICO'
catalog_amico = 'test'

outpath = "/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/amico/" + catalog_amico + "/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

#richness and mass cuts
min_richness = 0
inpath = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/amico/map_detections_refined_noBuffer_all.fits'

amico_data = Table.read(inpath) 
#print(cat_amico)
print(amico_data.colnames)

####General catalog properties
print("Number of AMICO clusters = ", len(amico_data))

#create clevar catalog
#c1 = ClCatalog('Cat1', ra=wazp_data['ra'], dec=wazp_data['dec'], z=wazp_data['redshift'], mass = wazp_data['NGALS'], id=wazp_data['ID'])
c1 = Table([amico_data['ID'],amico_data['Xphys'],amico_data['Yphys'],amico_data['Zphys'],amico_data['LAMB']],names=('id','ra','dec','z','mass'))

#add members
#c1.add_members(id=wazp_members_data['ID_g'], id_cluster=wazp_members_data['ID_CLUSTER'], ra=wazp_members_data['ra'], dec=wazp_members_data['dec'], pmem=wazp_members_data['PMEM'])
#c1_members = Table([wazp_members_data['ID_g'],wazp_members_data['ID_CLUSTER'],wazp_members_data['ra'],wazp_members_data['dec'],wazp_members_data['redshift'],wazp_members_data['PMEM']],names=('id','id_cluster','ra','dec','z','pmem'))
#print(c1)
#print(c1.members)
#c1.write(outpath + 'ClCatalog.fits', overwrite=True)
#c1.members.write(outpath + 'ClCatalog_members.fits', overwrite=True)
print(c1)
#print(c1_members)
c1.write(outpath + 'Catalog.fits', overwrite=True)
#c1_members.write(outpath + 'Catalog_members.fits', overwrite=True)

sys.exit()
