#!/usr/bin/env python
# coding: utf-8

###import
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import sys
import os
import shutil

inpath = "/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/"
inpath_2 = "/sps/lsst/groups/clusters/cluster_comparison_project/cosmoDC2/wazp/6948/"
outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/debug/wazp/"

#if os.path.exists(outpath):
#     shutil.rmtree(outpath)
os.makedirs(outpath,exist_ok=True)
print('outpath = ' + outpath)

#c1 = Table.read(inpath+'cosmoDC2/cosmoDC2_v1.1.4_small/Catalog_members.fits')
##c1 = Table.read(inpath+'wazp/6980/Catalog.fits')
#c2 = Table.read(inpath+'wazp/6980/Catalog_members.fits')
#c3 = Table.read(inpath_2+'Catalog.fits')
#c3 = Table.read(inpath+'redmapper/full/cosmoDC2_v1.1.4_redmapper_v0.8.1/Catalog.fits')
#c3 = Table.read(inpath+'amico/test/Catalog.fits')
c1 = Table.read('/sps/lsst/groups/clusters/cluster_comparison_project/debug/redmapper/Catalog.fits')
c2 = Table.read('/sps/lsst/groups/clusters/cluster_comparison_project/debug/wazp/Catalog.fits')
c3 = Table.read('/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/Catalog.fits')

print(c1.info)
print(c2.info)
print(c3.info)
#sys.exit()
#cuts
#c1=c1[c1['z']<1.15]
#c1 = c1[c1['snr']>4]
c1=c1[c1['mass']>10]
c2=c2[c2['mass']>10]
c3=c3[c3['mass']>10]
print(c2)
configuration = 'cosmoDC2 (full area): WaZP cluster redshift'
label_1 = 'redMaPPer'
label_2 = 'WaZP'
label_3 = 'AMICO'

#footprint check
#c3 = Table.read('/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/cosmoDC2/with_m200c/cosmoDC2_v1.1.4/Catalog_members.fits')
#c3 = Table.read('/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/cosmoDC2/m200c_gt_13.0/cosmoDC2_v1.1.4_image/Catalog_members.fits')
#c4 = Table.read('/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/cosmoDC2/with_m200c/cosmoDC2_v1.1.4/Catalog.fits')
#print(c3.info)
#print(c4.info)
#c3_cut_1=c3[c3['mag_i']>25]
#c3_cut_2=c3[c3['mag_i']<25]
#plt.plot(c4['ra'], c4['dec'], 'xr', label='halos')
#plt.scatter(c4['ra'],c4['dec'], marker='.',color = 'blue', s=0.02, alpha=0.3, label='halos')
#plt.xlim([45,80])
#plt.ylim([-50, -20])
#plt.xlabel('ra')
#plt.ylabel('dec')
#plt.title('cosmoDC2 halos (halo_mass>10**13)')
#plt.savefig(outpath+"cosmoDC2_halos.png", bbox_inches='tight')
#plt.close()

#try a nice density map
#fig, (ax1, ax2) = plt.subplots(1, 2)
#xbins = np.linspace(45, 80, 100)
#ybins = np.linspace(-50, -20, 100)
#fig, ax1 = plt.subplots()
#w, h = 2819, 100
#d_scale = [[0 for x in range(w)] for y in range(h)]
#d_scale = [20 for x in range(w)]
#for i in range(0,100):
#    d_scale[i]=100
#    for j in range(0,100):
#        #print('----' + str(i))
#        #print(j)
#        d_scale[i][j]=100
#print(len(xbins))
#print(d_scale)
#heatmap, xedges, yedges = np.histogram2d(c4['ra'], c4['dec'], bins=(xbins,ybins), weights=d_scale)
#heatmap, xedges, yedges = np.histogram2d(c4['ra'], c4['dec'], bins=100)
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
#ax2.set_xlim([45, 80])
#ax2.set_ylim([-50, -20])
#ax2.set_ylabel("DEC")
#ax2.set_xlabel("RA")
#fig.savefig(outpath+"map_clusters_smoothed.png", bbox_inches='tight')
#sys.exit()

#mass
#bin_range = [0,60]
#nbins = 60
#plt.figure()
#plt.hist(c1['mass'], range=bin_range, bins=nbins, label=label_1, histtype='step', color = 'black')
#plt.hist(c2['mass'], range=bin_range, bins=nbins, label=label_2, histtype='step', color = 'red')
#plt.hist(c3['mass'], range=bin_range, bins=nbins, label=label_3, histtype='step', color = 'blue')
#plt.xlabel("alg. richness");
#plt.ylabel("clusters")
#plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
#plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
#plt.legend(title = '', loc='upper right')
#plt.title(configuration)
#plt.savefig(outpath+'mass.png', bbox_inches='tight')
#plt.close() 

#redshift
bin_range = [0,2.0]
nbins = 20
plt.figure()
plt.hist(c1['z'], range=bin_range, bins=nbins, label=label_1, histtype='step', color = 'red', density=True)
plt.hist(c2['z'], range=bin_range, bins=nbins, label=label_2, histtype='step', color = 'black', density=True)
plt.hist(c3['z'], range=bin_range, bins=nbins, label=label_3, histtype='step', color = 'blue', density=True)
#plt.hist(c3_cut_2['z'], range=bin_range, bins=nbins, label='galaxies mag_i<25', histtype='step', color = 'red', density=True)
#plt.hist(c3_cut_2['z'], range=bin_range, bins=nbins, label='galaxies mag_i>25', histtype='step', color = 'black', density=True)
plt.xlabel("redshift");
plt.ylabel("% / 0.05 dz")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.legend(title = '', loc='upper right')
plt.title('cosmoDC2')
plt.savefig(outpath+'redshift.png', bbox_inches='tight')
plt.close() 
#sys.exit()

#cluster density versus richness
sky_area_sq_deg = 440
mag_bins = np.linspace(10, 100, 90)
cdf1 = np.searchsorted(c1["mass"], mag_bins, sorter=c1["mass"].argsort())
cdf1 = len(c1)-cdf1
cdf2 = np.searchsorted(c2["mass"], mag_bins, sorter=c2["mass"].argsort())
cdf2 = len(c2)-cdf2
cdf1 = np.searchsorted(c1["mass"], mag_bins, sorter=c1["mass"].argsort())
cdf1 = len(c1)-cdf1
cdf3 = np.searchsorted(c3["mass"], mag_bins, sorter=c3["mass"].argsort())
cdf3 = len(c3)-cdf3
g1, = plt.semilogy(mag_bins, cdf1 / sky_area_sq_deg, color = 'red')
g2, = plt.semilogy(mag_bins, cdf2 / sky_area_sq_deg, color = 'black')
g3, = plt.semilogy(mag_bins, cdf3 / sky_area_sq_deg, color = 'blue')
plt.xlabel("alg. richness");
plt.ylabel("Density (cl/deg2)");
plt.legend([g1,g2,g3],[label_1, label_2, label_3], title = 'cosmoDC2', loc='upper right')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.5', color='grey')
plt.savefig(outpath+"cluster_density.png", bbox_inches='tight')
plt.close()
sys.exit()

###################debug
c1=c1[c1['mass']>25]
c2=c2[c2['mass']>20]
d1=len(c1)/440
d2=len(c2)/440
print(c1)
print(c2)
print(d1)
print(d2)
sys.exit()
