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
from astropy.modeling import models, fitting

###cluster_validation
#from cluster_validation.opening_catalogs_functions import *
#from cluster_validation.wazp_functions import *

###lsst
#import lsst.afw.geom as afw_geom

###check noise.fits
#file='/sps/lsst/groups/clusters/amico_validation_project/catalogs/130323/cosmoDC2_small_photoz_flexzboost_v0/noise.fits'
#myfits=fits.open(file)
#myfits.info()
#print(myfits[0].header)
#imarr = myfits[0].data
#print(imarr)

#image_data = fits.getdata(file, ext=0)
#plt.figure(figsize=[5,10])
#plt.imshow(image_data, vmin=0, vmax=100000, cmap='hot', origin='lower')
#plt.colorbar()
#plt.savefig("photoz.png", bbox_inches='tight')
#plt.close()

#cosmoDC2
if(False):
    #file='/sps/lsst/groups/clusters/amico_validation_project/catalogs/cosmoDC2_photoz_flexzboost/v0/9554/galaxies.fits'
    file_1='/sps/lsst/users/tguillem/web/clusters/catalogs/cosmoDC2_photoz_flexzboost/v1/9554/galaxies.fits'
    galaxy_data_1 = Table.read(file_1)
    galaxy_data_1 = galaxy_data_1[galaxy_data_1['mag_i_photoz']>25]
    #galaxy_data_1 = galaxy_data_1[galaxy_data_1['mag_i_photoz']<25]
    file_2='/sps/lsst/users/tguillem/web/clusters/catalogs/cosmoDC2_photoz_flexzboost/v1/9425/galaxies.fits'
    galaxy_data_2 = Table.read(file_2)
    #['galaxy_id','ra','dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y', 'mag_g_photoz', 'mag_i_photoz', 'mag_r_photoz','mag_z_photoz', 'mag_y_photoz','photoz_mode', 'photoz_mean', 'photoz_mask', 'photoz_pdf', 'photoz_median', 'photoz_odds', 'photoz_mode_ml_red_chi2', 'photoz_mode_ml']

if(True):
    file_1='/sps/lsst/users/tguillem/web/clusters/catalogs/DC2_small_photoz_flexzboost/v0/4027/galaxies.fits'
    galaxy_data_1 = Table.read(file_1)
    galaxy_data_1 = galaxy_data_1[galaxy_data_1['redshift']>0.001]
    galaxy_data_2 = galaxy_data_1[galaxy_data_1['mag_i']<25.2]
    galaxy_data_3 = galaxy_data_1[galaxy_data_1['mag_y']<24.3]
    #galaxy_data_1 = galaxy_data_1[galaxy_data_1['redshift']>1.1]
    #galaxy_data_1 = galaxy_data_1[galaxy_data_1['redshift']>1.3]
    
fig = plt.figure(figsize=(12,8))
plt.hist(galaxy_data_1['redshift'], density=False, range=(0,2), bins = 400, label="zmode", histtype='step', color = 'black');
plt.hist(galaxy_data_2['redshift'], density=False, range=(0,2), bins = 400, label="zmode (mag_i<25.2)", histtype='step', color = 'red');
plt.hist(galaxy_data_3['redshift'], density=False, range=(0,2), bins = 400, label="zmode (mag_y<24)", histtype='step', color = 'blue');
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.xlabel("z");
plt.ylabel("galaxies");
#plt.legend([g1,g2], ["dr3", "cosmoDC2"], title = 'Photo-z', loc='upper right')
plt.legend(title = 'DC2 (4027) FlexZBoost', loc='upper right')
outpath = "/sps/lsst/users/tguillem/web/clusters/debug/"
plt.savefig(outpath+"redshift_gt_25.png", bbox_inches='tight')
plt.close()

#mag checks
#read mstar values
t_param=Table.read('/sps/lsst/users/tguillem/web/clusters/mstar/plots/m_star.fits')


plt.figure()
plt.scatter(galaxy_data_1['mag_y'], galaxy_data_1['redshift'], color = 'blue', s=0.01, alpha=0.3, label='galaxies')
plt.xlim([15, 30])
plt.ylim([0, 3])
plt.ylabel('redshift')
plt.xlabel('mag_i')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'scatter.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.scatter(galaxy_data_1['redshift'], galaxy_data_1['mag_y'], color = 'blue', s=0.01, alpha=0.3, label='galaxies')
plt.ylim([18, 28])
plt.xlim([0, 1.8])
plt.xlabel('redshift')
plt.ylabel('mag_y')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'scatter_mag_redshift.png', bbox_inches='tight')
plt.close()

#mag-mag plot
plt.figure()
plt.scatter(galaxy_data_1['mag_i'], galaxy_data_1['mag_y'], color = 'blue', s=0.01, alpha=0.3, label='galaxies')
plt.ylim([18, 28])
plt.xlim([18, 28])
plt.xlabel('mag_i')
plt.ylabel('mag_y')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'scatter_mag_y_vs_mag_i.png', bbox_inches='tight')
plt.close()

sys.exit()

###GCR catalogs import
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter
from GCRCatalogs.helpers.base_filters import sample_filter, partition_filter 
from GCRCatalogs import GCRQuery

###Access object table with GCRCatalogs
#GCRCatalogs.get_public_catalog_names()
#print(GCRCatalogs.get_available_catalog_names(include_default_only=False))
#obj_cat = GCRCatalogs.load_catalog("dc2_object_run2.2i_dr6_v2")

###Object Table schema
#print('obj_cat quantities')
#print(sorted(obj_cat.list_all_quantities()))

###Truth
#extragalactic_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small')
#galaxy_data_all = Table(extragalactic_cat.get_quantities(['galaxy_id', 'ra', 'dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y'])[0:10000])
#print('1000000 galaxies selected')
ra_min, ra_max = 64, 65
dec_min, dec_max = -40, -39
query = GCRQuery(
        'ra >= {}'.format(ra_min),
        'ra < {}'.format(ra_max),
        'dec >= {}'.format(dec_min),
        'dec < {}'.format(dec_max),
    )
#query = GCRCatalogs.GCRQuery('('ra > {}'.format(ra_min),'ra < {}'.format(ra_max),'dec > {}'.format(dec_min),'dec < {}'.format(dec_max))
#galaxy_data_all = Table(extragalactic_cat.get_quantities(['galaxy_id', 'ra', 'dec', 'halo_mass', 'halo_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y'],[query]))
#print(galaxy_data_all)
#extragalactic_cat = GCRCatalogs.load_catalog('skysim5000_v1.1.1_small')
#extragalactic_cat = GCRCatalogs.load_catalog('dc2_truth_run2.2i_summary_tract_partition')
#print('extragalactic_cat quantities')
#print(sorted(extragalactic_cat.list_all_quantities()))

###Truth catalog
#truth_cat = GCRCatalogs.load_catalog('dc2_truth_run2.2i_star_truth_summary')
#truth_cat = GCRCatalogs.load_catalog('dc2_truth_run2.2i_summary_tract_partition')
#print('truth_cat quantities')
#print(sorted(truth_cat.list_all_quantities()))
###Dump variable definitions
#for qty in ['cosmodc2_hp', 'cosmodc2_id', 'dec', 'flux_g', 'flux_g_noMW', 'flux_i', 'flux_i_noMW', 'flux_r', 'flux_r_noMW', 'flux_u', 'flux_u_noMW', 'flux_y', 'flux_y_noMW', 'flux_z', 'flux_z_noMW', 'host_galaxy', 'id', 'is_pointsource', 'is_variable', 'mag_g', 'mag_g_noMW', 'mag_i', 'mag_i_noMW', 'mag_r', 'mag_r_noMW', 'mag_u', 'mag_u_noMW', 'mag_y', 'mag_y_noMW', 'mag_z', 'mag_z_noMW', 'patch', 'ra', 'redshift', 'tract', 'truth_type']:
#    info_dict = truth_cat.get_quantity_info(qty)
#    print(qty,info_dict)

###Truth-match catalog
#truth_match_cat = GCRCatalogs.load_catalog('dc2_object_run2.2i_dr6_v2_with_addons')
#print('truth_match_cat quantities')
#print(sorted(truth_match_cat.list_all_quantities()))

###Photo-z catalog
#photoz_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small_with_photozs_v1')
#print('photoz_cat quantities')
#print(sorted(q for q in photoz_cat.list_all_quantities()if q.startswith('photoz_')))
#print(sorted(q for q in photoz_cat.list_all_quantities()))

###Define truth cuts
#cosmoDC2_cuts = [
    #is_galaxy,
    #bright,
#    GCRQuery('mag_i<24')  # Select objects that have i-band cmodel magnitudes
#]

###Truth selection
#columns_truth = ["ra", "dec", "mag_i", "redshift"]
#truth = extragalactic_cat.get_quantities(
#    quantities=columns_truth,
#    filters=cosmoDC2_cuts
#)

###Data selection
#is_extended = GCRQuery('extendedness == 1')  # Extended objects (primarily galaxies)
#clean = GCRQuery('clean')  # The source has no flagged pixels (interpolated, saturated, edge, clipped...) 
                           # and was not skipped by the deblender
#true_galaxy = GCRQuery('truth_type == 1')
#true_star = GCRQuery('truth_type == 2')
#true_SN = [GCRQuery('truth_type !=1'),GCRQuery('truth_type !=2')]
#is_good_match = GCRQuery('is_good_match == 1')
#true_galaxy_cosmoDC2 = GCRQuery('truth_id == 1')

#galaxy_cuts = [
#    clean,#is_good_match,
#    is_extended, 
    #true_galaxy,
    #GCRQuery((np.isfinite, 'mag_i_cModel')),  # Select objects that have i-band cmodel magnitudes
#    GCRQuery('mag_i_cModel<24')
#]

#star_cuts = [
#    clean,is_good_match,
#    #is_extended, 
#    true_star,
#    GCRQuery((np.isfinite, 'mag_i_cModel')),  # Select objects that have i-band cmodel magnitudes
#]

#is_extended_cuts = [
#    clean,is_good_match,
    #is_extended, 
    #true_SN,
#    GCRQuery('truth_type!=1'),
    #GCRQuery('truth_type !=2'),
#    GCRQuery((np.isfinite, 'mag_i_cModel')),  # Select objects that have i-band cmodel magnitudes
#]

###Quantities to store
#columns = ["ra", "dec", "photoz_mean", 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y']

###Photo-z catalog
photoz_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small_with_photozs_v1')#, config_overwrite={'healpix_pixels': ['9559',  '9686',  '9687', '9814']})
#photoz_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small_with_photozs_flexzboost_v1')
#photoz_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image_with_photozs_v1')
#photoz_cat = GCRCatalogs.load_catalog('dc2_object_run2.2i_dr6_with_addons', config_overwrite={'tracts': [3638]})
print('photoz_cat quantities')
print(sorted(q for q in photoz_cat.list_all_quantities()if q.startswith('photoz_')))
#print(sorted(q for q in photoz_cat.list_all_quantities()))

query = GCRQuery(
            'ra >= {}'.format(ra_min),
            'ra < {}'.format(ra_max),
            'dec >= {}'.format(dec_min),
            'dec < {}'.format(dec_max),
            #'mag_i<24',
            )
#data = Table(photoz_cat.get_quantities(columns,[query]))

#data = photoz_cat.get_quantities(['photoz_mean','photoz_mask', 'mag_i', 'mag_z'], native_filters=['healpix_pixel == 9559'])
data = photoz_cat.get_quantities(['redshift', 'photoz_mode', 'photoz_mean','photoz_mask', 'mag_i', 'mag_z'], native_filters=partition_filter('healpix_pixel',['9559',  '9686',  '9687', '9814']))
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

#,[query])
#galaxy_data_all = Table(photoz_cat.get_quantities('mag_i'))
#print(galaxy_data_all)
#galaxy_data_all = Table(photoz_cat.get_quantities('photoz_mean'))
#photoz_mask=galaxy_data_all['photoz_mask']
#galaxy_data_all = galaxy_data_all[photoz_mask]
#galaxy_data_all = Table(photoz_cat.get_quantities(['ra', 'dec', 'photoz_mean', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y'],[query]))
#galaxy_data_all.rename_column('photoz_mean','redshift')

#    native_filters=tract_filter([3830])#, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
#)

#data_galaxy_cut = photoz_cat.get_quantities(
#    quantities=columns, 
#    filters=galaxy_cuts 
#    #native_filters=tract_filter([3830])#, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
#)

#extragalactic_cut = extragalactic_cat.get_quantities(
#    quantities=columns_truth,
#    filters=cosmoDC2_cuts,
#    ) 

#Data_star_cut = Truth_match_cat.Get_quantities(
#    Quantities=Columns, 
#    Filters=Star_cuts, 
#    Native_filters=Tract_filter([3830, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
#)

#Data_extended_cut = Truth_match_cat.Get_quantities(
#    Quantities=Columns, 
#    Filters=Is_extended_cuts, 
#    Native_filters=Tract_filter([3830, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
#)

#sys.exit()

###mag selection
#print('********Start mag selection********')
#data_brightness_i = ascii.read("istar.asc")
#data_brightness_z = ascii.read("zstar.asc")
#filter_arr = []
#filter_arr2 = []
#for element in galaxy_data_all:
#    index_i = mstar_i(element['redshift'])
#    index_z = mstar_z(element['redshift'])
#    if ( (element['mag_i'] < data_brightness_i[index_i][1]+2) and (element['mag_i']<25.47) ):
##        #if ( element['mag_i']<25.47 ):
#        filter_arr.append(True)
#    else:
#        filter_arr.append(False)
#        
#    if ( (element['mag_z'] < data_brightness_z[index_z][1]+2) and (element['mag_z']<24.46) ):
#            #if ( element['mag_z']<24.46 ):
#            filter_arr2.append(True)
#    else:
#        filter_arr2.append(False)
#galaxy_data_1 = galaxy_data_all[filter_arr]
#galaxy_data_2 = galaxy_data_all[filter_arr2]

galaxy_data_1 = galaxy_data_all[galaxy_data_all['mag_i']<25.5]
galaxy_data_2 = galaxy_data_all[galaxy_data_all['mag_z']<24.5]

#SAVE


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
outpath = "/sps/lsst/users/tguillem/DESC/desc_may_2021/desc-data-portal/notebooks/dc2/plots/photoz/bpz/"
#sys.mkdir(outpath);
plt.savefig(outpath+"photoz_mode.png", bbox_inches='tight')
plt.close()
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
