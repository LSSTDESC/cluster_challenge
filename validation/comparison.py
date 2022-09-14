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
outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/debug/cosmoDC2/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

c1 = Table.read(inpath+'wazp/6980/Catalog.fits')
c2 = Table.read(inpath+'redmapper/full/cosmoDC2_v1.1.4_redmapper_v0.8.1/Catalog.fits')
#c3 = Table.read(inpath+'amico/test/Catalog.fits')
print(c1.info)
print(c2.info)
#cuts
c1=c1[c1['z']<1.15]
#c2=c2[c2['mass']>20]

configuration = 'cosmoDC2'
label_1 = 'WaZP'
label_2 = 'redMaPPer'
label_3 = 'AMICO'

#mass
bin_range = [0,60]
nbins = 60
plt.figure()
plt.hist(c1['mass'], range=bin_range, bins=nbins, label=label_1, histtype='step', color = 'black')
plt.hist(c2['mass'], range=bin_range, bins=nbins, label=label_2, histtype='step', color = 'red')
#plt.hist(c3['mass'], range=bin_range, bins=nbins, label=label_3, histtype='step', color = 'blue')
plt.xlabel("alg. richness");
plt.ylabel("clusters")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.legend(title = '', loc='upper right')
plt.title(configuration)
plt.savefig(outpath+'mass.png', bbox_inches='tight')
plt.close() 

#redshift
bin_range = [0,1.6]
nbins = 32
plt.figure()
plt.hist(c1['z'], range=bin_range, bins=nbins, label=label_1, histtype='step', color = 'black')
plt.hist(c2['z'], range=bin_range, bins=nbins, label=label_2, histtype='step', color = 'red')
#plt.hist(c3['z'], range=bin_range, bins=nbins, label=label_3, histtype='step', color = 'blue')
plt.xlabel("redshift");
plt.ylabel("clusters / 0.05 dz")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.legend(title = '', loc='upper right')
plt.title(configuration)
plt.savefig(outpath+'redshift.png', bbox_inches='tight')
plt.close() 

#cluster density versus richness
sky_area_sq_deg = 440
mag_bins = np.linspace(5, 100, 95)
cdf1 = np.searchsorted(c1["mass"], mag_bins, sorter=c1["mass"].argsort())
cdf1 = len(c1)-cdf1
cdf2 = np.searchsorted(c2["mass"], mag_bins, sorter=c2["mass"].argsort())
cdf2 = len(c2)-cdf2
g1, = plt.semilogy(mag_bins, cdf1 / sky_area_sq_deg, color = 'black')
g2, = plt.semilogy(mag_bins, cdf2 / sky_area_sq_deg, color = 'red')
plt.xlabel("alg. richness");
plt.ylabel("Density (cl/deg2)");
#plt.legend([g1,g2,g3,g4],["dr6", "dr6 (truth_type=1)", "dr6 (truth_type=2)", "dr6 (ext=1)"], title = 'DC2 (9 tracts) / is_good_match', loc='upper left')
#plt.legend([g1,g2,g3],["dr6", "dr6 (truth_type=1)", "dr6 (truth_type!=1)"], title = 'DC2 (9 tracts) / is_good_match', loc='upper left')
plt.legend([g1,g2],[label_1, label_2], title = 'cosmoDC2', loc='upper right')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.5', color='grey')
plt.savefig(outpath+"cluster_density.png", bbox_inches='tight')
plt.close()

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
