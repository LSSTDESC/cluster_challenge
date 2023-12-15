#!/usr/bin/env python
# coding: utf-8
print('Matching Initialisation') 
###import
import numpy as np
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import sys
import os

###clevar
import clevar
from clevar.catalog import ClCatalog
from clevar.match import ProximityMatch
from clevar.match import MembershipMatch
from clevar.match import get_matched_pairs
from clevar.match import output_matched_catalog
from clevar.cosmology import AstroPyCosmology
from clevar.match import output_catalog_with_matching


#---Parameters selection---#
##---Matching type(s)---##
param = str(sys.argv[2])
print('Using matching :', param)

##---Select catalogs to match---##
catalog = str(sys.argv[1])
cut = str(sys.argv[3])
#cut = ''
print('Opening catalogs : ', catalog, cut)



#---Paths---#

inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/"

#cosmoDC2
#amico_path_i = inpath + 'AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/'
#amico_path_y = inpath + 'AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_y/'
#amico_path = amico_path_i
#DC2
amico_path = inpath + "AMICO/amico_cats/DC2/mag_y/"

rm_path = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/redmapper/full_pmem/cosmoDC2_v1.1.4_redmapper_v0.8.1/'
#rm_path = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/RedMapper/' #Uncomment when cut on redMaPPer is performed

#cdc2_path = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/cosmoDC2/m200c_gt_13.0/cosmoDC2_v1.1.4_image/'
#cdc2_path = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/CosmoDC2/cosmoDC2_photoz_flexzboost/m14/" #Uncomment when cut on cosmoDC2 is performed

dc2_path = ""

outpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/"


#---Reading catalogs---#

if catalog == 'amico_cosmoDC2':
    print('AMICO CosmoDC2 matching')
    path1 = amico_path + 'cuts/'
    path2 = cdc2_path
    cut1 = cut
    cut2 = ''
    c1 = ClCatalog.read(amico_path + 'cuts/Catalog' + cut1 + '.fits', 'c1', full = True)
    c1.read_members(amico_path + 'cuts/Catalog_members' + cut1 + '.fits', full = True)
    c2 = ClCatalog.read(dc2_path + 'Catalog' + cut2 + '.fits', 'c2', full = True)
    c2.read_members(dc2_path + 'Catalog_members' + cut2 + '.fits', full = True)
    matching = 'amico_cosmoDC2/mag_i/w_rich/m14/cuts/'
    
elif catalog == 'amico_redmapper':
    print('AMICO RedMapper matching')
    path1 = amico_path
    path2 = rm_path
    cut1 = '_sn5'
    cut2 = '_ls12'
    c1 = ClCatalog.read(amico_path + 'Catalog' + cut1 + '.fits', 'c1', full = True)
    c1.read_members(amico_path + 'Catalog_members' + cut1 + '.fits', full = True)
    c2 = ClCatalog.read(rm_path + 'Catalog' + cut2 + '.fits', 'c2', full = True)
    c2.read_members(rm_path + 'Catalog_members' + cut2 + '.fits', full = True)  
    matching = 'amico_redmapper/mag_i/'
    
elif catalog == 'redmapper_cosmoDC2':
    path1 = rm_path
    path2 = dc2_path
    cut1 = ''
    cut2 = ''
    print('RedMapper CosmoDC2 matching')
    c1  = ClCatalog.read(rm_path + 'Catalog' + cut1 + '.fits', 'c2', full = True)
    c1.read_members(rm_path + 'Catalog_members' + cut1 + '.fits', full = True)
    c2 = ClCatalog.read(dc2_path + 'Catalog' + cut2 + '.fits', 'c2', full = True)
    c2.read_members(dc2_path + 'Catalog_members' + cut2 + '.fits', full = True)
    matching = 'redmapper_cosmoDC2/'
    
elif catalog == 'amico_amico':
    print('AMICO AMICO match')
    c1 = ClCatalog.read(amico_path_i + 'Catalog.fits', 'c1', full = True)
    c1.read_members(amico_path_i + 'Catalog_members.fits', full = True)
    c2 = ClCatalog.read(amico_path_y + 'Catalog.fits', 'c2', full = True)
    c2.read_members(amico_path_y + 'Catalog_members.fits', full = True)
    matching = 'amico_amico/'

else:
    print('Catalog selection is wrong.')
    sys.exit()

    
#---Condition---#
"""
#c2 = c2[c2['z']<1.15]
#c2.members = c2.members[c2.members['z']<1.15]
"""


#---Create final output folder---#

outpath += matching + ''
if not os.path.exists(outpath):
    os.makedirs(outpath)
print('outpath = ', outpath)


#---Matchings---#


##---Proximity matching---##

if param == 'p_matching':
    
    print('Proximity matching is processing')
    
    #Define catalogs without members
    c1_raw = c1.raw()
    c2_raw = c2.raw()
    
    mt = ProximityMatch()
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
    #mt.match_from_config(c1_raw, c2_raw, match_config, cosmo=cosmo)
    mt.match_from_config(c1_raw, c2_raw, match_config, cosmo=cosmo)
    c1_raw.write(outpath + 'c1_p' + cut1 + cut2 + '.fits', overwrite=True)
    c2_raw.write(outpath + 'c2_p' + cut1 + cut2 + '.fits', overwrite=True)
    c1_m = c1_raw[c1_raw['mt_cross'] != None]
    c2_m = c2_raw[c2_raw['mt_cross'] != None]
    c1_m.read_members(path1 + 'Catalog_members' + cut1 + '.fits', full = True)
    c2_m.read_members(path2 + 'Catalog_members' + cut2 + '.fits', full = True)  
    c1_m.members.write(outpath + 'c1_p_members' + cut1 + cut2 + '.fits', overwrite=True)
    c2_m.members.write(outpath + 'c2_p_members' + cut1 + cut2 + '.fits', overwrite=True)
    output_matched_catalog(outpath + 'c1_p' + cut1 + cut2 + '.fits', outpath + 'c2_p' + cut1 + cut2 + '.fits', outpath + 'output_catalog_p' + cut1 + cut2 + '.fits', c1_raw, c2_raw, matching_type='cross', overwrite = True)

##---Membership matching---##

elif param == 'mb_matching' :
    
    print('Membership matching is performing')
    
    fshare = 0.00
    print('fshare = ', fshare)
    outpath += str(fshare) + '/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    mt = MembershipMatch()
    match_config = {
      'type': 'cross', # options are cross, cat1, cat2
      'preference': 'shared_member_fraction', # other options are more_massive, angular_proximity or redshift_proximity
      'minimum_share_fraction1' : fshare,
      'minimum_share_fraction2' : fshare,
      'match_members_kwargs': {'method':'id'},
      }
    mt.match_from_config(c1, c2, match_config)
    c1.write(outpath + 'c1_mb.fits', overwrite=True)
    c2.write(outpath + 'c2_mb.fits', overwrite=True)
    mt.save_matches(c1, c2, out_dir=outpath, overwrite=True) 
    output_matched_catalog(outpath + 'c1_mb.fits', outpath + 'c2_mb.fits', outpath + 'output_catalog_mb.fits', c1, c2, matching_type='cross', overwrite = True)

    
sys.exit()