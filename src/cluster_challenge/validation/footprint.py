#!/usr/bin/env python
# coding: utf-8

###import
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import sys
import os
import shutil
import glob

inpath = "/sps/lsst/users/tguillem/web/clusters/catalogs/DC2_photoz_flexzboost/v0/"
#c1 = Table.read(inpath+'galaxies.fits')
#print(c1)
file_to_read=inpath+"*/galaxies.fits"
file_list=glob.glob(file_to_read)
#print(file_list)

outpath='/sps/lsst/users/tguillem/web/clusters/debug/'
plt.figure()
#for ifile in range(3):
for ifile in range(len(file_list)):
    
    c1 = Table.read(file_list[ifile])
    print(file_list[ifile] + ' : ' + str(len(c1)))
    #print(c1)

    plt.scatter(c1['ra'],c1['dec'], marker='.',color = 'blue', s=0.002, alpha=0.3)

plt.xlim([45,80])
plt.ylim([-50, -20])
plt.xlabel('ra')
plt.ylabel('dec')
plt.title('DC2 galaxies')
plt.savefig(outpath+"DC2.png", bbox_inches='tight')
plt.close()
