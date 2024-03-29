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
from clevar import ClCatalog
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

outpath = "/sps/lsst/groups/clusters/cluster_comparison_project/debug/"

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
c1 = Table([cluster_data['cluster_id'],cluster_data['ra'],cluster_data['dec'],cluster_data['redshift'],cluster_data['richness'],cluster_data['redshift_true_cg']],names=('id','ra','dec','z','mass','z_true'))

#add log(mass) column
c1.add_column(1.0, name='log_mass', index=6)
for i in range(0,len(c1)):
     c1['log_mass'][i]=np.log10(c1['mass'][i])

#add members
#c1.add_members(id=member_data['id_member'], id_cluster=member_data['cluster_id_member'], ra=member_data['ra_member'], dec=member_data['dec_member'])
c1_members = Table([member_data['id_member'],member_data['cluster_id_member'],member_data['ra_member'],member_data['dec_member'],member_data['mag_r_lsst_member'],member_data['mag_i_lsst_member'],member_data['mag_z_lsst_member']],names=('id','id_cluster','ra','dec','mag_r', 'mag_i', 'mag_z'))
#add correct pmem
pmem = member_data["p_member"] * member_data["pfree_member"] * member_data["theta_i_member"] * member_data["theta_r_member"] 
c1_members.add_column(pmem, name='pmem', index=4)

#restrict members to halos above richness requirement
filter_arr = []
for element in c1_members:
     if (element['id_cluster'] in c1['id']):
          filter_arr.append(True)
     else:
          filter_arr.append(False)
c1_members = c1_members[filter_arr]

#print(c1)
#print(c1.members)
#c1.write(outpath + 'ClCatalog.fits', overwrite=True)
#c1.members.write(outpath + 'ClCatalog_members.fits', overwrite=True)
print(c1)
print(c1_members)
c1.write(outpath + 'Catalog.fits', overwrite=True)
c1_members.write(outpath + 'Catalog_members.fits', overwrite=True)
sys.exit()

#try a nice density map
w, h = 2819, 100
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
heatmap, xedges, yedges = np.histogram2d(cluster_data['ra'], cluster_data['dec'], bins=(xbins,ybins), weights=d_scale)
##heatmap, xedges, yedges = np.histogram2d(cluster_data['ra_cen_0'], cluster_data['dec_cen_0'], bins=100)
##ax1.imshow(heatmap, interpolation='none', cmap='jet')
#im = ax1.pcolormesh(xbins, ybins, heatmap.T, cmap='jet')
#fig.colorbar(im, ax=ax1)
#fig.savefig(outpath+"map_clusters.png", bbox_inches='tight')
fig, ax2 = plt.subplots()
from astropy.convolution.kernels import Gaussian2DKernel
from astropy.convolution import convolve
##im = ax2.imshow(convolve(heatmap, Gaussian2DKernel(x_stddev=3,y_stddev=3)), interpolation='none', cmap='jet')
im = ax2.pcolormesh(xbins, ybins, convolve(heatmap, Gaussian2DKernel(x_stddev=2,y_stddev=2)).T, cmap='jet')
fig.colorbar(im, ax=ax2)
ax2.set_xlim([74, 49])
ax2.set_ylim([-46, -26])
ax2.set_ylabel("Dec")
ax2.set_xlabel("RA")
fig.savefig(outpath+"map_clusters_smoothed.png", bbox_inches='tight')

sys.exit()
