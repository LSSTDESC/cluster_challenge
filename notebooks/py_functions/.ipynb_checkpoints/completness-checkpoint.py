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
from scipy.optimize import minimize

###clevar
import clevar
from clevar.catalog import ClCatalog
from clevar.match import ProximityMatch
from clevar.match import get_matched_pairs
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances
from clevar.match_metrics.recovery import ClCatalogFuncs as r_cf
from clevar.match import output_matched_catalog

#matching_folder = '/sps/lsst/groups/clusters/redmapper_validation_project/cosmoDC2_v1.1.4/extragal/after_matching/v0/'
matching_folder = '/sps/lsst/groups/clusters/wazp_validation_project/cosmoDC2_v1.1.4/extragal/after_matching/v0/'

outpath_base = '/pbs/home/t/tguillem/web/clusters/cluster_challenge/selection_function/completeness_fit/'

###completeness

def completness_func1D(c1, c2, bin1, bin2, x_param): #give binnning from your params ex : bin1 = np.linspace(.2,1.2,9) & bin2 = [10**13,10**14,10**15]
    if x_param == 'z':
        param1, param2, param3, param4 = 'cat1_redshift', 'redshift', 'cat1_halo_mass', 'halo_mass'
    elif x_param == 'mass':
        param1, param2, param3, param4 = 'cat1_halo_mass', 'halo_mass', 'cat1_redshift', 'redshift'
    bin_range = [min(bin1), max(bin1)]
    nbins_x = len(bin1)-1
    compl = np.empty([len(bin2),nbins_x])
    for i in range(0,len(bin2)-1):
        print(i)
        cut1 = bin2[i]
        cut2 = bin2[i+1]
        filter1 = np.logical_and(c1.data[param3] > cut1, c1.data[param3] < cut2)
        filter2 = np.logical_and(c2.data[param4] > cut1, c2.data[param4] < cut2)
        c_halos_matched = c1[filter1]
        c_halos = c2.data[filter2]
        h_r_halos_matched = np.histogram(c_halos_matched[param1], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
        h_r_halos  = np.histogram(c_halos[param2], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
        compl[i] = np.divide(h_r_halos_matched[0], h_r_halos[0], where=(h_r_halos[0]!=0))
    return compl
    

def completness_func2D(c1, c2, x_bins, y_bins):
    #2D : z-M plane
    nbins_x = len(x_bins)
    nbins_y = len(y_bins)
    xbin_range = [min(x_bins), max(x_bins)]
    ybin_range = [min(y_bins), max(y_bins)]
    h2_z_halos_matched = np.histogram2d(c1['cat1_redshift'], c1['cat1_halo_mass'], bins=(x_bins, y_bins), 
                                        range=[xbin_range, ybin_range], normed=None, weights=None, density=None)
    h2_z_halos = np.histogram2d(c2.data['redshift'], c2.data['halo_mass'], bins=(x_bins, y_bins),
                                range=[xbin_range, ybin_range], normed=None, weights=None, density=None)
    number_of_match = h2_z_halos_matched[0]
    number_of_halo = h2_z_halos[0]
    compl_2d = (number_of_match/number_of_halo)
    return compl_2d
    
print('Completness module charged w 2d func')