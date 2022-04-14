#!/usr/bin/env python
# coding: utf-8

###import
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

###clevar
import clevar
from clevar.catalog import ClCatalog
from clevar.match import ProximityMatch
from clevar.match import MembershipMatch
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances
from clevar.cosmology import AstroPyCosmology
from clevar.match import output_catalog_with_matching

inpath = "/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/"
outpath = "/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/after_matching/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

#select the catalogs to match
wazp_cosmoDC2 = False
redmapper_cosmoDC2 = True
wazp_redmapper = False

if wazp_cosmoDC2 == True:
     #c1 = ClCatalog.read_full(inpath+'wazp/6685/ClCatalog.fits')
     #c1_members = ClCatalog.read_full(inpath+'wazp/6685/ClCatalog_members.fits')
     #c2 = ClCatalog.read_full(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/ClCatalog.fits')
     #c2_members = ClCatalog.read_full(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/ClCatalog_members.fits')
     c1  = ClCatalog.read(inpath+'wazp/6685/Catalog.fits', 'c1', id='id', ra='ra', dec='dec', z='z', mass='mass')
     #c1_members = ClCatalog.read(inpath+'wazp/6685/Catalog_members.fits', 'c1_members', id='id', id_cluster='id_cluster', ra='ra', dec='dec', z='z', pmem='pmem')
     c1.read_members(inpath+'wazp/6685/Catalog_members.fits',id='id', id_cluster='id_cluster', ra='ra', dec='dec', z='z', pmem='pmem')
     c2  = ClCatalog.read(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog.fits', 'c2', id='id', ra='ra', dec='dec', z='z', mass='mass')
     #c2_members = ClCatalog.read(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog_members.fits', 'c2_members', id='id', id_cluster='id_cluster', ra='ra', dec='dec', z='z', pmem='pmem')
     c2.read_members(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog_members.fits', id='id', id_cluster='id_cluster', ra='ra', dec='dec', z='z', pmem='pmem')
elif redmapper_cosmoDC2 == True:
     c1  = ClCatalog.read(inpath+'redmapper/cosmoDC2_v1.1.4_redmapper_v0.8.1/Catalog.fits', 'c1', id='id', ra='ra', dec='dec', z='z', mass='mass')
     #c1_members = ClCatalog.read(inpath+'redmapper/cosmoDC2_v1.1.4_redmapper_v0.8.1/Catalog_members.fits', 'c1_members', id='id', id_cluster='id_cluster', ra='ra', dec='dec',pmem='pmem')
     c1.read_members(inpath+'redmapper/cosmoDC2_v1.1.4_redmapper_v0.8.1/Catalog_members.fits', id='id', id_cluster='id_cluster', ra='ra', dec='dec',pmem='pmem')
     c2  = ClCatalog.read(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog.fits', 'c2', id='id', ra='ra', dec='dec', z='z', mass='mass')
     #c2_members = ClCatalog.read(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog_members.fits', 'c2_members', id='id', id_cluster='id_cluster', ra='ra', dec='dec', z='z', pmem='pmem')
     c2.read_members(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog_members.fits', id='id', id_cluster='id_cluster', ra='ra', dec='dec', z='z', pmem='pmem')
else:
     print('Catalog selection is wrong.')
     sys.exit()

#define catalogs without members
c1_raw = c1.raw()
c2_raw = c2.raw()

###perform proximity matching
if True:
     mt = ProximityMatch()

     #old way but working     
     #mt_config1 = {'delta_z':0.05,
     #              'match_radius': '1 mpc',
     #              'cosmo':AstroPyCosmology()}
     #mt_config2 = {'delta_z':0.05,
     #              'match_radius': '1 mpc',
     #              'cosmo':AstroPyCosmology()}
     #mt.prep_cat_for_match(c1, **mt_config1)
     #mt.prep_cat_for_match(c2, **mt_config2)
     #mt.multiple(c1, c2)
     #mt.multiple(c2, c1)
     #mt.unique(c1, c2, preference='angular_proximity')
     #mt.unique(c2, c1, preference='angular_proximity')
     #c1.cross_match()
     #c2.cross_match()

     match_config = {
          'type': 'cross', # options are cross, cat1, cat2
          'which_radius': 'max', # Case of radius to be used, can be: cat1, cat2, min, max
          'preference': 'angular_proximity', # options are more_massive, angular_proximity or redshift_proximity
          'catalog1': {'delta_z':.05,
                       'match_radius': '1 mpc'
                       },
          'catalog2': {'delta_z':.05,
                       'match_radius': '1 mpc'
                       }
          }
     cosmo = AstroPyCosmology()
          
     mt.match_from_config(c1_raw, c2_raw, match_config, cosmo=cosmo)
     c1_raw.cross_match()
     c2_raw.cross_match()
     c1_raw.write(outpath + 'c1.fits', overwrite=True)
     c2_raw.write(outpath + 'c2.fits', overwrite=True)
     #mt.save_matches(c1_raw, c2_raw, out_dir=outpath, overwrite=True)
     
     #to print summary
     #mt.load_matches(c1, c2, out_dir=outpath)

###perform member matching
if False:
     mt = MembershipMatch()
     match_config = {
          'type': 'cross', # options are cross, cat1, cat2
          'preference': 'shared_member_fraction', # other options are more_massive, angular_proximity or redshift_proximity
          'minimum_share_fraction': 0,
          'match_members_kwargs': {'method':'id'},
          }
     mt.match_from_config(c1, c2, match_config)
     c1.write(outpath + 'c1.fits', overwrite=True)
     c2.write(outpath + 'c2.fits', overwrite=True)

sys.exit()
