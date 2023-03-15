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

def purity_func1D(c1, c2, bin1, bin2, bin_param): #give binnning from your params ex : bin1 = np.linspace(.2,1.2,9) & bin2 = [10**13,10**14,10**15]
    if bin_param == 'z':
        param1, param2, param3, param4 = 'cat2_redshift', 'redshift', 'cat2_mass', 'mass' 
    elif bin_param == 'richness':
        param1, param2, param3, param4 = 'cat2_mass', 'mass', 'cat2_redshift', 'redshift'
    bin_range = [min(bin1), max(bin1)]
    nbins_x = len(bin1)-1
    purity = np.empty([len(bin2),nbins_x])
    for i in range(0,len(bin2)-1):
        print(i)
        cut1 = bin2[i]
        cut2 = bin2[i+1]
        filter1 = np.logical_and(c1.data[param1] > cut1, c1.data[param1] < cut2)
        c_clusters_matched = c1[filter1]
        filter2 = np.logical_and(c2.data[param2] > cut1, c2.data[param2] < cut2)
        c_clusters = c2.data[filter2]
        h_r_clusters_matched = np.histogram(np.log(c_clusters_matched[param3]), bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
        h_r_clusters  = np.histogram(np.log(c_clusters[param4]), bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
        purity[i] = np.divide(h_r_clusters_matched[0], h_r_clusters[0], where=(h_r_clusters[0]!=0))
    return purity

def purity_func2D(c1, c2, x_bins, y_bins):
    #2D : z-M plane
    nbins_x = len(x_bins)
    nbins_y = len(y_bins)
    xbin_range = [min(x_bins), max(x_bins)]
    ybin_range = [min(y_bins), max(y_bins)]
    h2_z_halos_matched = np.histogram2d(c1['cat1_z'], c1['richness'], bins=(x_bins, y_bins), 
                                        range=[xbin_range, ybin_range], normed=None, weights=None, density=None)
    h2_z_halos = np.histogram2d(c2.data['z'], c2.data['richness'], bins=(x_bins, y_bins),
                                range=[xbin_range, ybin_range], normed=None, weights=None, density=None)
    number_of_match = h2_z_halos_matched[0]
    number_of_halo = h2_z_halos[0]
    purity_2d = (number_of_match/number_of_halo)
    return purity_2d

print('purity module charged w 2d func')
