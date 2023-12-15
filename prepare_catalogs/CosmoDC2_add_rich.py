#!/usr/bin/env python
# coding: utf-8
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii

###cluster_validation
from amico_functions import *


cdc2_path = "/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/cosmoDC2/m200c_gt_13.0/cosmoDC2_v1.1.4_image/"
mstar_path = "/sps/lsst/users/tguillem/web/clusters/mstar/"
outpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/CosmoDC2/cosmoDC2_photoz_flexzboost/m13/"
#mass cuts
min_halo_mass = 10**13 #Msun
#max_halo_mass = 10**13.5
#gc and gc_truth: catalog objects / others are just tables
halo_data = Table.read(cdc2_path + 'Catalog.fits')
halo_data = halo_data[(halo_data['mass']>=min_halo_mass)]#*(halo_data['mass']<=max_halo_mass)]

#halo members
galaxy_data_all = Table.read(cdc2_path + 'Catalog_members.fits')

print("Number of elements in the truth catalog = ", len(galaxy_data_all))
print("Number of halos in the truth catalog = ", len(halo_data))
print(galaxy_data_all, halo_data)
###mag selection
print('********Start mag selection********')
data_brightness_i = ascii.read(mstar_path + "DC2_i_star.dat")
data_brightness_y = ascii.read(mstar_path + "DC2_y_star.dat")
filter_arr = []
filter_arr2 = []
for element in galaxy_data_all:
    index_i = mstar_i(element['z'])
    index_y = mstar_y(element['z'])
    if ( (element['mag_i'] < data_brightness_i[index_i][1]+2) and (element['mag_i']<25.47) ):
        filter_arr.append(True)
    else:
        filter_arr.append(False)
    if ( (element['mag_y'] < data_brightness_y[index_y][1]+2) and (element['mag_y']<25.47) ):
        filter_arr2.append(True)
    else:
        filter_arr2.append(False)

galaxy_data_1 = galaxy_data_all[filter_arr]
galaxy_data_2 = galaxy_data_all[filter_arr2]
print('********End mag selection********')

print('********Start richness computation********')
richness_i = []
richness_y = []
for halo in halo_data:
    #_i
    members_i = galaxy_data_1[galaxy_data_1['id_cluster']==halo['id']]
    #print(members_i)
    richness_i.append(len(members_i))
    #_y
    members_y = galaxy_data_2[galaxy_data_2['id_cluster']==halo['id']]
    #print(members_z)
    richness_y.append(len(members_y))
     
halo_data.add_column(richness_i, index=5, name = 'NGALS_i')
halo_data.add_column(richness_y, index=5, name = 'NGALS_y')
print(halo_data)
print('********End richness computation********') 
halo_data.write(outpath + 'Catalog.fits')
galaxy_data_1.write(outpath+'Catalog_members.fits')
