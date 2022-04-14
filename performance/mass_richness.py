#!/usr/bin/env python
# coding: utf-8

###import
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
import numpy as np
from scipy.optimize import curve_fit
import sys
import os
import shutil
import pickle
import math

###clevar
import clevar
from clevar.catalog import ClCatalog
from clevar.match import ProximityMatch
from clevar.match import get_matched_pairs
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances

matching_folder = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/after_matching/'

##########select case
#catalog1 = 'wazp.fits'
#catalog1 = 'redmapper.fits'
#catalog2 = 'halos.fits'
#catalog2 = 'wazp.fits'
#catalog2 = 'redmapper.fits'
catalog1 = 'c1.fits'
catalog2 = 'c2.fits'
##########

outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/debug/"
if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)

#load c1 and c2
c1 = ClCatalog.read_full(matching_folder + catalog1)
c2 = ClCatalog.read_full(matching_folder + catalog2)
#print(c1.data)
#print(c2.data)

#restrict to matched pairs
#mt1, mt2 = get_matched_pairs(c1, c2, 'cross', None, None) 


#mass_density_metrics
plt.figure()
info = scaling.mass_density_metrics(
            c2, c1, 'cross', ax_rotation=45,
            add_fit=True, fit_bins1=8)
plt.savefig(outpath+'scaling_mass_density_metrics.png', bbox_inches='tight')
plt.close()

sys.exit()
