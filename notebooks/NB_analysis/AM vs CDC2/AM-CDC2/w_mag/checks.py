#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, hstack, vstack
from clevar.catalog import ClCatalog
from numpy.ma import masked
import sys

inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/amico_cosmoDC2/mag_i/w_rich/"
outpath = inpath
#am_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/Catalog_members.fits' 
#cdc_inpath = "/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/cosmoDC2/m200c_gt_13.0/cosmoDC2_v1.1.4_image/Catalog_members.fits"

print('reading tables')
c1 = ClCatalog.read(inpath + 'c1_p.fits', 'c1', full = True)
c1 = c1[c1['mt_cross']!= None]
c2 = ClCatalog.read(inpath + 'c2_p.fits', 'c2', full = True)
c2 = c2[c2['mt_cross']!= None]
c1.read_members(inpath + 'c1_p_members.fits', full = True)
c2.read_members(inpath + 'c2_p_members.fits', full = True)
#c2.members = c2.members[c2.members['mag_i']<=25.3]
print("Done, now it's just plotting", max(c1.members['mag_i']), max(c2.members['mag_i']), len(c1), len(c2))

print(len(c1.members), len(c2.members))
