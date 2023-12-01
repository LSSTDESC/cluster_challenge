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
inpath_members = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/amico_map_associations/10071_map_associations.fits'

amico_data = Table.read(inpath) 
print(amico_data.colnames)

amico_members_data = Table.read(inpath_members)
print(amico_members_data.colnames)
#amico_members.pprint_all()

####General catalog properties
print("Number of AMICO clusters = ", len(amico_data))
print("Number of AMICO members (several entries per member) = ", len(amico_members_data))

#create clevar catalog
#c1 = ClCatalog('Cat1', ra=wazp_data['ra'], dec=wazp_data['dec'], z=wazp_data['redshift'], mass = wazp_data['NGALS'], id=wazp_data['ID'])
c1 = Table([amico_data['uid'],amico_data['Xphys'],amico_data['Yphys'],amico_data['Zphys'],amico_data['LAMBSTAR']],names=('id','ra','dec','z','mass'))
#fix by hand the id issue
#c1.add_column(1.0, name='id', index=1)
#for i in range(0,len(c1)):
#     c1['id'][i]=i+1
#c1['id']=c1['id'].astype(int)

#add members
c1_members = Table([amico_members_data['mb_id'],amico_members_data['halo_id'], amico_members_data['prob']], names=('id','id_cluster','pmem'))
c1_members = c1_members[c1_members['pmem']>0.5]
print(c1)
print(c1_members)
c1.write(outpath + 'Catalog.fits', overwrite=True)
c1_members.write(outpath + 'Catalog_members.fits', overwrite=True)

sys.exit()
