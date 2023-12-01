#!/usr/bin/env python
# coding: utf-8

# Modified from:
# Rubin LSST DESC DC2: Accessing Object Table with GCRCatalogs

###Import necessary packages
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii

###GCR catalogs import
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter
from GCRCatalogs.helpers.base_filters import sample_filter, partition_filter 
from GCRCatalogs import GCRQuery

sys.path.append('../prepare_catalogs/')
from wazp_functions import *

###Photo-z catalog
photoz_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small_with_photozs_v1')#, config_overwrite={'healpix_pixels': ['9559',  '9686',  '9687', '9814']})
#photoz_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small_with_photozs_flexzboost_v1')
#photoz_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image_with_photozs_v1')
#photoz_cat = GCRCatalogs.load_catalog('dc2_object_run2.2i_dr6_with_addons', config_overwrite={'tracts': [3638]})
print('photoz_cat quantities')
print(sorted(q for q in photoz_cat.list_all_quantities()if q.startswith('photoz_')))
#print(sorted(q for q in photoz_cat.list_all_quantities()))

#data = photoz_cat.get_quantities(['photoz_mean','photoz_mask', 'mag_i', 'mag_z'], native_filters=['healpix_pixel == 9559'])
data = photoz_cat.get_quantities(['redshift', 'photoz_mode', 'photoz_mean','photoz_mask', 'mag_i', 'mag_z'], native_filters=partition_filter('healpix_pixel',['9559']))
                                                                                                                                                              #,  '9686',  '9687', '9814']))
# | healpix_pixel == 9686'])#,  9686,  9687,  9814]'])
#native_filters={'healpix_pixel': [9559,  9686,  9687,  9814]})# OR 'healpix_pixel==9687'])
#config_overwrite={'healpix_pixels': [8786, 8787, 8788]})
###HACK because of the size issue
redshift=data['redshift']
photoz_mean=data['photoz_mean']
photoz_mode=data['photoz_mode']
photoz_mask=data['photoz_mask']
mag_i=data['mag_i']
mag_z=data['mag_z']
#photoz_mean=photoz_mean[photoz_mask]
redshift=redshift[photoz_mask]
mag_i=mag_i[photoz_mask]
mag_z=mag_z[photoz_mask]
galaxy_data_all = Table([redshift,photoz_mean,photoz_mode,mag_i,mag_z],names=('redshift_true', 'photoz_mean', 'photoz_mode', 'mag_i', 'mag_z'))
#galaxy_data_all = galaxy_data_all[0:1000000]
print(galaxy_data_all)

###mag selection
#print('********Start mag selection********')
#data_brightness_i = ascii.read("istar.asc")
data_brightness_z = ascii.read("/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/prepare_catalogs/mstar_files/zstar.asc")
filter_arr = []
for element in galaxy_data_all:
         index_z = mstar_z(element['redshift_true'])
         if ( (element['mag_z'] < data_brightness_z[index_z][1]+2) and (element['mag_z']<24.46)):
             filter_arr.append(True)
         else:
             filter_arr.append(False)
galaxy_data_all = galaxy_data_all[filter_arr]

galaxy_data_1 = galaxy_data_all[galaxy_data_all['mag_i']<30]
galaxy_data_2 = galaxy_data_all[galaxy_data_all['mag_z']<30]

###Photoz plot
#fig = plt.figure(figsize=(12,8))
#plt.hist(x, range = (0, 3), bins = 30)
#plt.hist(data['photoz_mean'], density=True, range = (0, 3), bins = 60, label="photo-z mode", histtype='step', color = 'purple');
#plt.hist(galaxy_data_all['redshift'], density=True,range = (0, 3), bins = 30, label="true z", histtype='step', color = 'black');
#plt.hist(extragalactic_cat['redshift'], density=True, range = (0, 3), bins = 30, label="true z", histtype='step', color = 'black');
#plt.hist(galaxy_data_2['redshift'], density=True, range = (0, 3), bins = 30, label="true z (m* z-band)", histtype='step', color = 'blue');
#plt.hist(galaxy_data_1['redshift'], density=True, range = (0, 3), bins = 30, label="true z (m* i-band)", histtype='step', color = 'red');
#plt.plot(cat.photoz_pdf_bin_centers, sumpdf*3.,label="summed $p(z)$",lw=2,c='r');
###HACK
#plt.hist(galaxy_data_all['redshift_true'], density=False,range = (0, 2), bins = 40, label="redshift", histtype='step', color = 'black');
#plt.hist(galaxy_data_all['redshift'], density=False,range = (0, 2), bins = 40, label="pz", histtype='step', color = 'purple');
plt.hist(galaxy_data_1['redshift_true'], density=False,range = (0, 2), bins = 40, label="redshift (mag_i<25.5)", histtype='step', color = 'black');
plt.hist(galaxy_data_1['photoz_mean'], density=False, range = (0, 2), bins = 40, label="pz_mean (mag_i<25.5)", histtype='step', color = 'red');
plt.hist(galaxy_data_1['photoz_mode'], density=False, range = (0, 2), bins = 40, label="pz_mode (mag_i<25.5)", histtype='step', color = 'blue');
#plt.hist(galaxy_data_2['redshift'], density=False, range = (0, 2), bins = 40, label="pz (mag_z<24.5)", histtype='step', color = 'blue');
#plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
# Customize the minor grid
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.5', color='grey')
plt.xlabel("z");
plt.ylabel("galaxies / 0.05 dz");
#plt.legend([g1,g2], ["dr3", "cosmoDC2"], title = 'Photo-z', loc='upper right') 
plt.legend(title = 'cosmoDC2_small (BPZ)', loc='upper right') 
outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/redshift/wazp_cosmoDC2/"
#sys.mkdir(outpath);
plt.savefig(outpath+"photoz_values.png", bbox_inches='tight')
plt.close()

plt.figure()
plt.scatter(galaxy_data_1['photoz_mode'],galaxy_data_1['photoz_mode']-galaxy_data_1['redshift_true'], marker='.',color = 'blue', s=0.5, alpha=0.3, label='galaxies')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0, 1.5])
plt.ylim([-0.3, 0.3])
plt.xlabel('zmode')
plt.ylabel('zmode-z_true')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'redshift_correlation.png', bbox_inches='tight')

sys.exit()











###Compute sky area
d_ra = data["ra"].max() - data["ra"].min()
d_dec = data["dec"].max() - data["dec"].min()
cos_dec = np.cos(np.deg2rad(np.median(data["dec"])))
sky_area_sq_deg = (d_ra * cos_dec) * d_dec
print(sky_area_sq_deg)
print('********Plot********')
###Data plot
mag_bins = np.linspace(14, 24, 200)
cdf = np.searchsorted(data["mag_i_cModel"], mag_bins, sorter=data["mag_i_cModel"].argsort())
cdf_galaxy_cut = np.searchsorted(data_galaxy_cut["mag_i_cModel"], mag_bins, sorter=data_galaxy_cut["mag_i_cModel"].argsort())
cdf_star_cut = np.searchsorted(data_star_cut["mag_i_cModel"], mag_bins, sorter=data_star_cut["mag_i_cModel"].argsort())
cdf_extended_cut = np.searchsorted(data_extended_cut["mag_i_cModel"], mag_bins, sorter=data_extended_cut["mag_i_cModel"].argsort())
g1, = plt.semilogy(mag_bins, cdf / sky_area_sq_deg, color = 'black')
g2, = plt.semilogy(mag_bins, cdf_galaxy_cut / sky_area_sq_deg, color = 'red')
#g3, = plt.semilogy(mag_bins, cdf_star_cut / sky_area_sq_deg, color = 'blue')
#g4, = plt.semilogy(mag_bins, cdf_extended_cut / sky_area_sq_deg, color = 'orange')
###cosmoDC2
Cdf = np.searchsorted(truth["mag_i"], mag_bins, sorter=truth["mag_i"].argsort())
###Compute sky area cosmoDC2
d_ra = truth["ra"].max() - truth["ra"].min()
d_dec = truth["dec"].max() - truth["dec"].min()
cos_dec = np.cos(np.deg2rad(np.median(truth["dec"])))
sky_area_sq_deg_cosmoDC2 = (d_ra * cos_dec) * d_dec
print(sky_area_sq_deg_cosmoDC2)
g5, = plt.semilogy(mag_bins, Cdf / sky_area_sq_deg_cosmoDC2, color = 'blue')
plt.xlabel("$i$-band mag");
plt.ylabel("Galaxies /deg2 /0.05mag");
#plt.legend([g1,g2,g3,g4],["dr6", "dr6 (truth_type=1)", "dr6 (truth_type=2)", "dr6 (ext=1)"], title = 'DC2 (9 tracts) / is_good_match', loc='upper left')
#plt.legend([g1,g2,g3],["dr6", "dr6 (truth_type=1)", "dr6 (truth_type!=1)"], title = 'DC2 (9 tracts) / is_good_match', loc='upper left')
plt.legend([g1,g2,g5],["dr6", "dr6 (ext=1)", "cosmoDC2"], title = 'DC2 (9 tracts)', loc='upper left')
#plt.grid();
#plt.grid(True, which="both", ls="-", )
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
# Customize the minor grid
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.5', color='grey')
###Define the output path for figures
outpath = "/sps/lsst/users/tguillem/DESC/config/tests_DC2/desc-data-portal/notebooks/dc2/plots/matching/"
#sys.mkdir(outpath);
plt.savefig(outpath+"dr6.png", bbox_inches='tight') 
plt.close()
print('********Plot saved********')

###CosmoDC2 plot
# Cdf = Np.Searchsorted(Truth["Mag_i_cModel"], Mag_bins, Sorter=Truth["Mag_i_cModel"].Argsort())
# Sky_area_sq_deg = 50;
# Plt.Semilogy(Mag_bins, Cdf / Sky_area_sq_deg)
# Plt.Xlabel("$I$-Band Mag");
# Plt.Ylabel("Galaxies /Deg2 /0.05mag");
# Plt.Grid();  # Add Grid To Guide Our Eyes
# ###Define The Output Path For Figures
# Outpath = "/Sps/Lsst/Users/Tguillem/Desc/Config/Tests_dc2/Desc-Data-Portal/Notebooks/Dc2/Plots/"
# Plt.Savefig(Outpath+"Truth.Png", Bbox_inches='Tight')


###Filtering out the stars using the matched catalog
  ###to use the matched catalog to filter out the stars


###Redshift distributions
