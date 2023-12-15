import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, hstack, vstack
from clevar.catalog import ClCatalog
from numpy.ma import masked
import sys
from clevar.match import output_matched_catalog
inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/amico_cosmoDC2/mag_i/"
c_merged12 = ClCatalog.read(inpath + 'output_catalog_p.fits')
am_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/Catalog.fits' 
cdc_inpath = "/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/cosmoDC2/m200c_gt_13.0/cosmoDC2_v1.1.4_image/Catalog.fits"
c1 = ClCatalog.read(inpath + 'c1_p.fits', 'c1', full = True)
c1.read_members('/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/Catalog_members.fits', full = True)
print(c1)
print(c1.members, len(c1.members))