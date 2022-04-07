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
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances

inpath = "/sps/lsst/users/tguillem/DESC/desc_may_2021/cluster_challenge/clevar_catalogs/"
outpath = "/sps/lsst/users/tguillem/DESC/desc_may_2021/cluster_challenge/clevar_catalogs/after_matching/"

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
     c1_members = ClCatalog.read(inpath+'wazp/6685/Catalog_members.fits', 'c1_members', id='id', id_cluster='id_cluster', ra='ra', dec='dec', z='z', pmem='pmem')
     c2  = ClCatalog.read(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog.fits', 'c2', id='id', ra='ra', dec='dec', z='z', mass='mass')
     c2_members = ClCatalog.read(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog_members.fits', 'c2_members', id='id', id_cluster='id_cluster', ra='ra', dec='dec', z='z', pmem='pmem')
elif redmapper_cosmoDC2 == True:
     c1  = ClCatalog.read(inpath+'redmapper/cosmoDC2_v1.1.4_redmapper_v0.8.1/Catalog.fits', 'c1', id='id', ra='ra', dec='dec', z='z', mass='mass')
     c1_members = ClCatalog.read(inpath+'redmapper/cosmoDC2_v1.1.4_redmapper_v0.8.1/Catalog_members.fits', 'c1_members', id='id', id_cluster='id_cluster', ra='ra', dec='dec',pmem='pmem')
     c2  = ClCatalog.read(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog.fits', 'c2', id='id', ra='ra', dec='dec', z='z', mass='mass')
     c2_members = ClCatalog.read(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog_members.fits', 'c2_members', id='id', id_cluster='id_cluster', ra='ra', dec='dec', z='z', pmem='pmem')
else:
     print('Catalog selection is wrong.')
     sys.exit()

#perform proximity matching
from clevar.match import ProximityMatch
mt = ProximityMatch()
     
from clevar.cosmology import AstroPyCosmology
mt_config1 = {'delta_z':0.05,
              'match_radius': '1 mpc',
              'cosmo':AstroPyCosmology()}
mt_config2 = {'delta_z':0.05,
              'match_radius': '1 mpc',
              'cosmo':AstroPyCosmology()}
mt.prep_cat_for_match(c1, **mt_config1)
mt.prep_cat_for_match(c2, **mt_config2)

mt.multiple(c1, c2)
mt.multiple(c2, c1)

mt.unique(c1, c2, preference='angular_proximity')
mt.unique(c2, c1, preference='angular_proximity')

c1.cross_match()
c2.cross_match()

c1.write(outpath + 'c1.fits', overwrite=True)
c2.write(outpath + 'c2.fits', overwrite=True)
mt.save_matches(c1, c2, out_dir=outpath, overwrite=True)
               
#to print summary
mt.load_matches(c1, c2, out_dir=outpath)
sys.exit()








#if os.path.exists('matching/' + outpath_matching):
#     shutil.rmtree('matching/' + outpath_matching)
#os.makedirs('matching/' + outpath_matching)

###WaZP: different format
#cat_wazp = fits.open('catalogs/wazp_cluster.fits')
#cat_wazp_table = cat_wazp[0]
#cat_wazp_original = Table.read('catalogs/wazp_cluster.fits')
#cat_wazp = 
#print(cat_wazp.info)
#richness and mass cuts
min_richness = 0
wazp_cat_name = 'catalogs/' + catalog_wazp + '/wazp_cluster.fits'
wazp_members_cat_name = 'catalogs/' + catalog_wazp + '/wazp_membership.fits'
#truth matching
#matching = np.load('catalogs/wazp_truth_matching.npy')
#print(matching)
#test DES
#wazp_cat_name = 'catalogs/y1a1_dnf_wazp_v5.0.11.5239+47_clusters.fits'
#wazp_members_cat_name = 'catalogs/y1a1_dnf_wazp_v5.0.11.5239+47_members.fits'
wazp_data, wazp_members_data = wazp_cat_open(wazp_cat_name, wazp_members_cat_name, min_richness)
####General catalog properties
print("Number of WaZP clusters = ", len(wazp_data))
print("Number of WaZP cluster members = ", len(wazp_members_data))

###Redmapper
####function to open truth and detection catalogs
#RM_cat_name = 'cosmoDC2_v1.1.4_redmapper_v0.8.0'
RM_cat_name = 'cosmoDC2_v1.1.4_redmapper_v0.8.1'
###catalog from images
#RM_cat_name = 'dc2_redmapper_run2.2i_dr6_wfd_v0.8.1'
#richness and mass cuts
min_richness = 0
cluster_data_all, member_data_all, gc = RM_cat_open(RM_cat_name, min_richness, cluster_only=True)
print("Number of Redmapper clusters = ", len(cluster_data_all))
print("Number of Redmapper cluster members = ", len(member_data_all))
print("Redmapper sky area = ", gc.sky_area, "deg2")
#restrict to small cosmoDC2 footprint
healpix_pixels = [9559,  9686,  9687,  9814,  9815,  9816,  9942,  9943, 10070, 10071, 10072, 10198, 10199, 10200, 10326, 10327, 10450]
nside = 32
filter_arr = []
all_healpix_pixels = []
for element in cluster_data_all:
     pix = hp.ang2pix(nside, element['ra'], element['dec'], lonlat=True)
     if ( pix in healpix_pixels ):
          filter_arr.append(True)
     else:
          filter_arr.append(False)
     if ( pix not in all_healpix_pixels ):
          all_healpix_pixels.append(pix)
cluster_data = cluster_data_all[filter_arr]
#cluster_data = cluster_data_all
#print(all_healpix_pixels)
#print(len(all_healpix_pixels))
#print(hp.pixelfunc.pix2ang(nside, 10070, nest=False, lonlat=True))
#print(hp.pixelfunc.max_pixrad(nside,degrees=True))
#print(hp.pixelfunc.get_interp_val(10070, 61.5, -40.5, nest=False, lonlat=True))
#hp_32 = HEALPix(nside=32)#, order='nested')
#index, dx, dy = hp_32.lonlat_to_healpix(10070,return_offsets=True)
#print(index)
#print(dx)
#print(dy)
#do ra/dec scans
ra_center=[]
dec_center=[]
ra_min=[]
dec_min=[]
ra_max=[]
dec_max=[]
#healpix_pixels = [10070]
print('****ra/dec of healpix pixels by hand****')
all_healpix_pixels=[]
for element in all_healpix_pixels:
     #get center
     v_ra_center,v_dec_center = hp.pixelfunc.pix2ang(nside, element, nest=False, lonlat=True) 
     v_ra_center=round(v_ra_center,3)
     v_dec_center=round(v_dec_center,3)
     #get extrema
     v_ra_min=v_ra_center
     v_ra_max=v_ra_center
     v_dec_min=v_dec_center
     v_dec_max=v_dec_center
     for index in range(0,4000):
          ra=v_ra_center-2+index*0.001
          dec=v_dec_center-2+index*0.001
          pix = hp.ang2pix(nside, ra, v_dec_center, lonlat=True)
          if(pix==element and ra<v_ra_min):
               v_ra_min=ra
          if(pix==element and ra>v_ra_max):
               v_ra_max=ra
          pix = hp.ang2pix(nside, v_ra_center, dec, lonlat=True)     
          if(pix==element and dec<v_dec_min):
               v_dec_min=dec
          if(pix==element and dec>v_dec_max):
               v_dec_max=dec     
     #append
     ra_center.append(v_ra_center)
     dec_center.append(v_dec_center)
     ra_min.append(round(v_ra_min,3))
     dec_min.append(round(v_dec_min,3))
     ra_max.append(round(v_ra_max,3))
     dec_max.append(round(v_dec_max,3))
print(healpix_pixels)
print(ra_center)
print(dec_center)
print(ra_min)
print(dec_min)
print(ra_max)
print(dec_max)
#Create a table
#healpix_table = Table([all_healpix_pixels, ra_center, dec_center, ra_min, dec_min, ra_max, dec_max],names=('healpix_pixel', 'ra_center', 'dec_center', 'ra_min', 'dec_min', 'ra_max', 'dec_max'))
#healpix_table.write('/sps/lsst/users/tguillem/DESC/desc_may_2021/desc-data-portal/notebooks/dc2/tables/map_cosmoDC2.fits')
#print(healpix_table)
#healpix_table.pprint_all() 
#sys.exit()

#not applied
#cluster_data = cluster_data_all
filter_arr = []
for element in member_data_all:
     pix = hp.ang2pix(nside, element['ra_member'], element['dec_member'], lonlat=True)
     if ( pix in healpix_pixels ):
          filter_arr.append(True)
     else:
          filter_arr.append(False)
member_data = member_data_all[filter_arr]
#not applied
#member_data = member_data_all

###check of the areas
###Basic visualization
#fig, ax = plt.subplots()
##fig.figure()
#ax.plot(cluster_data['ra_cen_0'],cluster_data['dec_cen_0'],'r.', markersize=2, label = 'RM clusters)')
#ax.plot(wazp_data['ra'],wazp_data['dec'],'b.', markersize=1, label = 'WaZP clusters')
#ax.set_xlim([75, 48])
#ax.set_ylim([-46, -24])
#ax.set_xlabel("ra")
#ax.set_ylabel("dec");
#ax.set_title("cosmoDC2 small")
#ax.locator_params(axis="x", nbins=27)
#ax.locator_params(axis="y", nbins=23)
#ax.grid(which='major', axis='both', linestyle='--', linewidth='1', color='grey')
##ax.add_patch(Rectangle((-44,72),20,14, edgecolor='black', facecolor='none'))
##box: ra 52-72 / dec -44 - -30
#pp1 = plt.Rectangle((52, -44),20,14,edgecolor='black',facecolor='none',lw=4)
#ax.add_patch(pp1)
##fig.legend(loc='best', framealpha=0.3)
#ax.legend(bbox_to_anchor = (1, 1), loc = 'best', prop = {'size': 15})
#fig.savefig(outpath+"clusters.png", bbox_inches='tight')
#fig.close()
#print('********Plot saved********')

#try a nice density map
#fig, (ax1, ax2) = plt.subplots(1, 2)
#xbins = np.linspace(49, 74, 100)
#ybins = np.linspace(-46, -26, 100)
#fig, ax1 = plt.subplots()
#w, h = 2819, 100
#d_scale = [[0 for x in range(w)] for y in range(h)] 
#d_scale = [20 for x in range(w)]
##for i in range(0,100):
##     d_scale[i]=100
#     #for j in range(0,100):
#          #print('----' + str(i))
#          #print(j)
#          #d_scale[i][j]=100
##print(len(xbins))
##print(d_scale)
#heatmap, xedges, yedges = np.histogram2d(cluster_data['ra_cen_0'], cluster_data['dec_cen_0'], bins=(xbins,ybins), weights=d_scale)
##heatmap, xedges, yedges = np.histogram2d(cluster_data['ra_cen_0'], cluster_data['dec_cen_0'], bins=100)
##ax1.imshow(heatmap, interpolation='none', cmap='jet')
#im = ax1.pcolormesh(xbins, ybins, heatmap.T, cmap='jet')
#fig.colorbar(im, ax=ax1)
#fig.savefig(outpath+"map_clusters.png", bbox_inches='tight')
#fig, ax2 = plt.subplots()
#from astropy.convolution.kernels import Gaussian2DKernel
#from astropy.convolution import convolve
##im = ax2.imshow(convolve(heatmap, Gaussian2DKernel(x_stddev=3,y_stddev=3)), interpolation='none', cmap='jet')
#im = ax2.pcolormesh(xbins, ybins, convolve(heatmap, Gaussian2DKernel(x_stddev=2,y_stddev=2)).T, cmap='jet')
#fig.colorbar(im, ax=ax2)
#ax2.set_xlim([74, 49])
#ax2.set_ylim([-46, -26]) 
#ax2.set_ylabel("Dec")
#ax2.set_xlabel("RA")
#fig.savefig(outpath+"map_clusters_smoothed.png", bbox_inches='tight')
##sys.exit()
#print("Number of Redmapper clusters small= ", len(cluster_data))
#sys.exit()
####cluster-halo matching with clevar
####perform matching versus cosmoDC2
##DC2_cat_name = 'cosmoDC2_v1.1.4'
DC2_cat_name = 'cosmoDC2_v1.1.4_small'
##richness and mass cuts
min_halo_mass = 10**13 #Msun
##gc and gc_truth: catalog objects / others are just tables
truth_data, gc_truth = DC2_cat_open(DC2_cat_name, min_halo_mass, cluster_only=False)

#halo table
halo_data = truth_data[truth_data['is_central']==True]

###read hdf5 files from skysim
##inpath = "/sps/lsst/users/ccombet/SkySim5000/hdf5/"
##with pd.HDFStore(os.path.join(inpath,f'skysim_halos_z=0-1.20_mfof_gt_1.00e+13_image.hdf5')) as store:
##     halo_data = store['skysim']
##     halo_metadata = store.get_storer('skysim').attrs.metadata
##     #print(halo_data['baseDC2/sod_halo_mass'])

###rename
##halo_data.rename(columns={'baseDC2/sod_halo_mass': 'M200c', 'richness': 'NGALS', 'richness_i': 'NGALS_i', 'richness_z': 'NGALS_z'}, inplace=True)
###fix M200c
##halo_data['M200c'] = halo_data['M200c']/0.71
##print(halo_data)
##halo_data=halo_data[halo_data['M200c']>10**13.5]
###convert to an astropy table
##halo_data=Table.from_pandas(halo_data)

#restrict to DC2 rectangle footprint: ra 52-72 / dec -44 - -30
##filter_arr = []
##ra_min, ra_max = 52, 72
##dec_min, dec_max = -44, -30
##for element in halo_data:
##    if(element['ra']>ra_min and element['ra']<ra_max and element['dec']>dec_min and element['dec']<dec_max):
##          filter_arr.append(True)
##     else: 
##          filter_arr.append(False)
#apply or not DC2 rectangle selection
##halo_data = halo_data[filter_arr]

#mask = truth_data['is_central']==True
#mask = np.logical_and(truth_data['is_central']==True,truth_data['halo_id']==127200143167)
#halo_data = truth_data[mask]
#restrict to DR6 footprint
#healpix_pixels = [9559,  9686,  9687,  9814,  9815,  9816,  9942,  9943, 10070, 10071, 10072, 10198, 10199, 10200, 10326, 10327, 10450]
#nside = 32
#filter_arr = []
#for element in halo_data_all:
#     pix = hp.ang2pix(nside, element['ra'], element['dec'], lonlat=True)
#     if ( pix in healpix_pixels ):
#          filter_arr.append(True)
#     else:
#          filter_arr.append(False)
#halo_data = halo_data_all[filter_arr]

#halo members
#galaxy_data_all = Table(gc_truth.get_quantities(['ra','dec', 'halo_mass', 'halo_id', 'galaxy_id', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y', 'baseDC2/is_on_red_sequence_gr', 'baseDC2/is_on_red_sequence_ri', 'is_central', 'hostHaloMass'],
##,#, 'baseDC2/sod_halo_mass'],
#                                                filters=['is_central==False', 'halo_mass > 10**13']))
# + ']')halo_mass > ' + str(min_halo_mass) '+']))# 'halo_id==127200143167']))

#print("Number of elements in the truth catalog = ", len(galaxy_data_all))
#print("Number of halos in the truth catalog = ", len(halo_data))
#print("Truth catalog sky area = ", gc_truth.sky_area, "deg2")
#print(halo_data)
#print(galaxy_data_all)

#restrict to z<1.15
#wazp_data = wazp_data[wazp_data['redshift']<1.15]
#cluster_data = cluster_data[cluster_data['redshift']<1.15]

###check redmapper unmatched clusters
#cluster_data_unmatched = Table.read('plots/matching/wazp_6685_rgt25_redmapper_rgt20_proximity/redmapper_unmatched.fits')
#list_id = []
#for element in cluster_data_unmatched:
#     if element['id'] not in list_id:
#          list_id.append(int(element['id']))
#print(list_id)
#restrict clusters to this list of id
#filter_arr = []
#for element in cluster_data:
#     if element['cluster_id'] in list_id:
#          print(element['cluster_id'])
#          filter_arr.append(True)
#     else:
#          filter_arr.append(False)
#cluster_data = cluster_data[filter_arr]
#print(cluster_data)

sys.exit()
#matching
from clevar.catalog import ClCatalog
if algo == 'WaZP':
     c1 = ClCatalog('Cat1', ra=wazp_data['ra'], dec=wazp_data['dec'], z=wazp_data['redshift'], mass = wazp_data['NGALS'], id=wazp_data['ID'])
elif algo == 'redMaPPer':
     c1 = ClCatalog('Cat1', ra=cluster_data['ra'], dec=cluster_data['dec'], z=cluster_data['redshift'], mass = cluster_data['richness'], id=cluster_data['cluster_id'])
else:
     print('Cluster algorithm not defined')
     sys.exit()

#c2 = ClCatalog('Cat2', ra=halo_data['ra'], dec=halo_data['dec'], z=halo_data['redshift'], mass=halo_data['halo_mass'], ngals_i = halo_data['NGALS_i'], ngals_z = halo_data['NGALS_z'])
c2 = ClCatalog('Cat2', ra=halo_data['ra'], dec=halo_data['dec'], z=halo_data['redshift'], mass=halo_data['halo_mass'], id=halo_data['halo_id'])
#c2 = ClCatalog('Cat2', ra=halo_data['ra'], dec=halo_data['dec'], z=halo_data['redshift'], mass=halo_data['M200c'], id=halo_data['halo_id'])
#WaZP-redMaPPer
#c1 = ClCatalog('Cat1', ra=wazp_data['ra'], dec=wazp_data['dec'], z=wazp_data['redshift'], mass = wazp_data['NGALS'], id=wazp_data['ID'])
#c2 = ClCatalog('Cat2', ra=cluster_data['ra'], dec=cluster_data['dec'], z=cluster_data['redshift'], mass = cluster_data['richness'], id=cluster_data['cluster_id'])
#print(c1)
#print(c2)

#create clevar member catalogs
#from clevar.catalog import MemCatalog
#m1 = MemCatalog('Mem1', id=wazp_members_data['ID_g'], id_cluster=wazp_members_data['ID_CLUSTER'], ra=wazp_members_data['ra'], dec=wazp_members_data['dec'], pmem=wazp_members_data['PMEM'])
##m2 = MemCatalog('Mem2', id=galaxy_data_all['galaxy_id'], id_cluster=galaxy_data_all['halo_id'], ra=galaxy_data_all['ra'], dec=galaxy_data_all['dec']) # pmem=input2_mem['PMEM'])
#m2 = MemCatalog('Mem2', id=member_data['id_member'], id_cluster=member_data['cluster_id_member'], ra=member_data['ra_member'], dec=member_data['dec_member']) # pmem=input2_mem['PMEM'])
#m1 = MemCatalog('Mem1', id=member_data['id_member'], id_cluster=member_data['cluster_id_member'], ra=member_data['ra_member'], dec=member_data['dec_member'])
#m2 = MemCatalog('Mem2', id=galaxy_data_all['galaxy_id'], id_cluster=galaxy_data_all['halo_id'], ra=galaxy_data_all['ra'], dec=galaxy_data_all['dec'])

#clean m2
#mask = [CL_ID in c2.id_dict for CL_ID in m2['id_cluster']]
#m2 = m2[mask]
#print(m1)
#print(m2)

#test
#c2 = ClCatalog('Cat2', id=input2['ID'], ra=input2['RA'], dec=input2['DEC'], z=input2['Z'], mass=input2['MASS'])
#print(c1)
#print(c2)
#save c1 and c2
#wazp
#c1.data['ra', 'dec', 'z', 'mass'].write('matching/' + outpath_matching + '/wazp.fits', overwrite=True)
#c2.data['ra', 'dec', 'z', 'mass'].write('matching/' + outpath_matching + '/halos.fits', overwrite=True)
#c1.data['ra', 'dec', 'z', 'mass'].write('matching/' + outpath_matching + '/redmapper.fits')
#redmapper
#c1.data['ra', 'dec', 'z', 'mass'].write('matching/' + outpath_matching + '/redmapper.fits', overwrite=True)
#c2.data['ra', 'dec', 'z', 'mass'].write('matching/' + outpath_matching + '/halos.fits', overwrite=True)
#sys.exit()

print('******************************')
#clean catalog
#mask = [CL_ID in c1.id_dict for CL_ID in m1['id_cluster']]
#m1 = m1[mask]
#print(m1)
#mask2 = [CL_ID in c1.id_dict for CL_ID in m2['id_cluster']]
#m2 = m2[mask2]
#print(m2)

#test member matching
#add members
#from clevar.catalog import MemCatalog
#c1.add_members(id=wazp_members_data['ID_g'], id_cluster=wazp_members_data['ID_CLUSTER'], ra=wazp_members_data['ra'], dec=wazp_members_data['dec'], pmem=wazp_members_data['PMEM'])
#c2.add_members(id=member_data['id_member'], id_cluster=member_data['cluster_id_member'], ra=member_data['ra_member'], dec=member_data['dec_member'])
#display(c1.members)
#display(c2.members)

#from clevar.match import MembershipMatch
#mt = MembershipMatch()

#step-by-step method
#mt.match_members(c1.members, c2.members, method='id')
#print('matching done: 1')
#print(mt.matched_mems)
#mt.fill_shared_members(c1, c2, m1, m2)
#print('matching done: 2')
#mt.multiple(c1, c2)
#mt.multiple(c2, c1)
#mt.unique(c1, c2, 'shared_member_fraction')
#mt.unique(c2, c1, 'shared_member_fraction')
#mt.cross_match(c1)
#mt.cross_match(c2)
#mt.save_matches(c1, c2, 'mem_pmem', overwrite=True)
#mt.load_matches(c1, c2, 'mem_pmem')
#print('member matching done')

#or with config object
#match_config = {
#     'type': 'cross', # options are cross, cat1, cat2
#     'preference': 'shared_member_fraction', # other options are more_massive, angular_proximity or redshift_proximity
#     'minimum_share_fraction': 0,
#     'match_members_kwargs': {'method':'id'},
#     }

#mt.match_from_config(c1, c2, m1, m2, match_config)
#mt.match_from_config(c1, c2, match_config)
#mt.save_matches(c1, c2, out_dir='matching/' + outpath_matching + '/', overwrite=True)
#to print summary
#mt.load_matches(c1, c2, out_dir='matching/' + outpath_matching + '/')
#mt.save_matches(c1, c2, 'mem_pmem', overwrite=True)
#mt.load_matches(c1, c2, 'mem_pmem')
##sys.exit()

if True:
     from clevar.match import ProximityMatch
     mt = ProximityMatch()
     
     from clevar.cosmology import AstroPyCosmology
     mt_config1 = {'delta_z':0.05,
                   'match_radius': '1 mpc',
                   'cosmo':AstroPyCosmology()}
     mt_config2 = {'delta_z':0.05,
                   'match_radius': '1 mpc',
                   'cosmo':AstroPyCosmology()}
     mt.prep_cat_for_match(c1, **mt_config1)
     mt.prep_cat_for_match(c2, **mt_config2)

     mt.multiple(c1, c2)
     mt.multiple(c2, c1)

     mt.unique(c1, c2, preference='angular_proximity')
     mt.unique(c2, c1, preference='angular_proximity')
     
     c1.cross_match()
     c2.cross_match()
     
if algo == 'WaZP':
     c1.write(outpath + 'wazp.fits', overwrite=True)
     c2.write(outpath + 'halos.fits', overwrite=True)
     #c2.write(outpath + 'redmapper.fits', overwrite=True)
if algo == 'redMaPPer': 
     c1.write(outpath + '/redmapper.fits', overwrite=True)
     c2.write(outpath + 'halos.fits', overwrite=True)
mt.save_matches(c1, c2, out_dir=outpath, overwrite=True)
               
#to print summary
mt.load_matches(c1, c2, out_dir=outpath)
#display(c1)
#display(c2)
#print(c1)
#print(c2)

sys.exit()
