#!/usr/bin/env python
# coding: utf-8

###import
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import sys
import os
import shutil

inpath = "/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/cosmoDC2/m200c_gt_13.0/cosmoDC2_v1.1.4/"
outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/debug/wazp/"

#if os.path.exists(outpath):
#     shutil.rmtree(outpath)
os.makedirs(outpath,exist_ok=True)
print('outpath = ' + outpath)

c1 = Table.read(inpath+'Catalog.fits')
print(c1.info)
#sys.exit()

plt.figure()
plt.scatter(c1['z'],c1['log_m200c'], marker='.',color = 'blue', s=1, alpha=0.3, label='halos')
plt.xlim([1.1,1.5])
plt.ylim([13, 15])
plt.xlabel('z')
plt.ylabel('log(m200c)')
plt.title('cosmoDC2 halos (m200c>10**13)')
plt.savefig(outpath+"cosmoDC2_halos.png", bbox_inches='tight')
plt.close()
sys.exit()

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
