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
algo = 'redMaPPer'
#catalog
#RM_cat_name = 'cosmoDC2_v1.1.4_redmapper_v0.8.0'
RM_cat_name = 'cosmoDC2_v1.1.4_redmapper_v0.8.1'
#catalog from images
#RM_cat_name = 'dc2_redmapper_run2.2i_dr6_wfd_v0.8.1'

outpath = "/sps/lsst/users/tguillem/DESC/desc_may_2021/cluster_challenge/clevar_catalogs/redmapper/" + RM_cat_name + "/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

#richness and mass cuts
min_richness = 0
cluster_data, member_data, gc = RM_cat_open(RM_cat_name, min_richness, cluster_only=True)
print("Number of Redmapper clusters = ", len(cluster_data))
print("Number of Redmapper cluster members = ", len(member_data))
print("Redmapper sky area = ", gc.sky_area, "deg2")

#restrict to small cosmoDC2 footprint
filter_cosmoDC2_small = False
if filter_cosmoDC2_small == True:
     healpix_pixels = [9559,  9686,  9687,  9814,  9815,  9816,  9942,  9943, 10070, 10071, 10072, 10198, 10199, 10200, 10326, 10327, 10450]
     nside = 32
     filter_arr = []
     for element in cluster_data:
          pix = hp.ang2pix(nside, element['ra'], element['dec'], lonlat=True)
          if ( pix in healpix_pixels ):
               filter_arr.append(True)
          else:
               filter_arr.append(False)
               if ( pix not in all_healpix_pixels ):
                    all_healpix_pixels.append(pix)
     cluster_data = cluster_data[filter_arr]

     filter_arr = []
     for element in member_data:
          pix = hp.ang2pix(nside, element['ra_member'], element['dec_member'], lonlat=True)
          if ( pix in healpix_pixels ):
               filter_arr.append(True)
          else:
               filter_arr.append(False)
     member_data = member_data[filter_arr]

#create clevar catalog
#c1 = ClCatalog('Cat1', ra=cluster_data['ra'], dec=cluster_data['dec'], z=cluster_data['redshift'], mass = cluster_data['richness'], id=cluster_data['cluster_id'])
c1 = Table([cluster_data['cluster_id'],cluster_data['ra'],cluster_data['dec'],cluster_data['redshift'],cluster_data['richness']],names=('id','ra','dec','z','mass'))

#add members
#c1.add_members(id=member_data['id_member'], id_cluster=member_data['cluster_id_member'], ra=member_data['ra_member'], dec=member_data['dec_member'])
c1_members = Table([member_data['id_member'],member_data['cluster_id_member'],member_data['ra_member'],member_data['dec_member']],names=('id','id_cluster','ra','dec'))
c1_members.add_column(1.0, name='pmem', index=4)

#print(c1)
#print(c1.members)
#c1.write(outpath + 'ClCatalog.fits', overwrite=True)
#c1.members.write(outpath + 'ClCatalog_members.fits', overwrite=True)
print(c1)
print(c1_members)
c1.write(outpath + 'Catalog.fits', overwrite=True)
c1_members.write(outpath + 'Catalog_members.fits', overwrite=True)

sys.exit()
