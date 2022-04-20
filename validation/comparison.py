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
outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/debug/amico/"

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)
print('outpath = ' + outpath)

c1 = Table.read(inpath+'wazp/6685/Catalog.fits')
c2 = Table.read(inpath+'redmapper/cosmoDC2_v1.1.4_redmapper_v0.8.1/Catalog.fits')
c3 = Table.read(inpath+'amico/test/Catalog.fits')
print(c1.info)
print(c2.info)

configuration = 'cosmoDC2 small'
label_1 = 'WaZP'
label_2 = 'redMaPPer'
label_3 = 'AMICO'

#mass
bin_range = [0,60]
nbins = 30
plt.figure()
plt.hist(c1['mass'], range=bin_range, bins=nbins, label=label_1, histtype='step', color = 'black')
plt.hist(c2['mass'], range=bin_range, bins=nbins, label=label_2, histtype='step', color = 'red')
plt.hist(c3['mass'], range=bin_range, bins=nbins, label=label_3, histtype='step', color = 'blue')
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
nbins = 16
plt.figure()
plt.hist(c1['z'], range=bin_range, bins=nbins, label=label_1, histtype='step', color = 'black')
plt.hist(c2['z'], range=bin_range, bins=nbins, label=label_2, histtype='step', color = 'red')
plt.hist(c3['z'], range=bin_range, bins=nbins, label=label_3, histtype='step', color = 'blue')
plt.xlabel("redshift");
plt.ylabel("clusters")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.legend(title = '', loc='upper right')
plt.title(configuration)
plt.savefig(outpath+'redshift.png', bbox_inches='tight')
plt.close() 

sys.exit()
