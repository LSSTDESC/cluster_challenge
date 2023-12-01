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
from scipy import optimize
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.integrate import simps
from numpy import exp, linspace, random
import emcee
import corner

###clevar
import clevar
from clevar.catalog import ClCatalog
from clevar.match import ProximityMatch
from clevar.match import get_matched_pairs
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances
from clevar.match_metrics.scaling import ClCatalogFuncs as s_cf
from clevar.match import output_matched_catalog

#matching_folder = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/batch/after_member_matching_fix/full/m200c/'
#matching_folder = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/after_matching/'
#/sps/lsst/groups/clusters/redmapper_validation_project/cosmoDC2_v1.1.4/extragal/after_matching/v0/
#'/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/batch/after_member_matching/m200c/'
#matching_folder = matching_folder + 'amico_cosmoDC2/'
#matching_folder = '/sps/lsst/groups/clusters/redmapper_validation_project/cosmoDC2_v1.1.4/extragal/after_matching/v0/'
matching_folder = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/after_matching/member_f0/wazp_cosmoDC2/'

##########select case
catalog1 = 'c1.fits'
catalog2 = 'c2.fits'
##########

outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/redshift/wazp_cosmoDC2/"
if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)

##########fit functions
def gauss(x, a, x0, sigma):
     return a*np.exp(-(x-x0)**2/(2*sigma**2))

def sum_2gauss(x, a1, x1, sigma1, a2, x2, sigma2):
     return gauss(x, a1, x1, sigma1)+gauss(x, a2, x2, sigma2)

def power(x,a,b):
     return a*np.power(x,b)

r0 = 35
z0 = 0.8
def mass_richness(m0, x1, x2, f, g):
     #return np.log10(m0) + f * np.log10(r/r0) + g * np.log10((1+z)/(1+z0))
     return np.log10(m0) + f * np.log10(x1/r0) + g * np.log10((1+x2)/(1+z0))

def mass_richness_1D(m0, x1, f):
     return m0 + f * np.log10(x1/r0)
##########

#load c1 and c2
c1 = ClCatalog.read_full(matching_folder + catalog1)
c2 = ClCatalog.read_full(matching_folder + catalog2)
print(c1.data[0])
print(c2.data[0])

#create a merged catalog for the cross-matched pairs
#output_matched_catalog(matching_folder+catalog1, matching_folder+catalog2,matching_folder+'output_catalog.fits', c1, c2, matching_type='cross', overwrite=True)
#c_merged = ClCatalog.read(matching_folder+'output_catalog.fits', 'merged',  z_cl='cat1_z', richness = 'cat1_mass', z_halo='cat2_z', mass = 'cat2_mass', log_mass = 'cat2_log_mass', m200c = 'cat2_m200c', log_m200c = 'cat2_log_m200c',)
#c_merged = ClCatalog.read(matching_folder+'output_catalog.fits', 'merged',  z='cat1_z', richness = 'cat1_mass', z_halo='cat2_z', mass = 'cat2_mass', log_mass = 'cat2_log_mass')
#c_merged = c_merged[c_merged.data['m200c']>10**14]
#c_merged = c_merged[c_merged.data['richness']>20]
#c_merged.add_column(1.0, name='log_richness', index=6)
#for i in range(0,len(c_merged)):
#          c_merged.data['log_mass'][i]=np.log10(c_merged.data['mass'][i])
#print(c_merged)
#c_merged.write(outpath + 'c_merged.fits', overwrite=True)
#sys.exit()

mt1, mt2 = get_matched_pairs(c1, c2, 'cross')
#correlation plots
plt.figure()
plt.scatter(mt1['z'],mt1['z']-mt2['z'], marker='.',color = 'blue', s=0.5, alpha=0.3, label='galaxies')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0, 1.5])
plt.ylim([-0.3, 0.3])
plt.xlabel('z_cl')
plt.ylabel('z_cl-z_halo')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'redshift_correlation.png', bbox_inches='tight')
#sys.exit()

#mass_density_metrics
plt.figure()
info = scaling.mass_density_metrics(
            c2, c1, 'cross', ax_rotation=45,
            add_fit=True, fit_bins1=6)
#print(info)
#info['ax'].set_xlabel('$halo_mass$')
#info['ax'].set_ylabel('LAMBSTAR')
#info['ax'].set_ylim(0,100)
#info['ax'].set_xlim(10**13.5,10**15)
#info['ax'].set_title(matching) 
plt.savefig(outpath+'scaling_mass_density_metrics.png', bbox_inches='tight')
plt.close()

#mass_density_metrics from merged catalog: does not work
#plt.figure()
#info = scaling.mass_density_metrics(
#            c_merged.data['cat2_mass'], c_merged.data['cat1_mass'], 'cross', ax_rotation=45,
#            add_fit=True, fit_bins1=8)
#plt.savefig(outpath+'scaling_mass_density_metrics.png', bbox_inches='tight')
#plt.close()

#log(richness) VS log(mass)
plt.figure()
#info = scaling.mass_density_metrics(
#            c2, c1, 'cross', ax_rotation=45,
#            add_fit=True, fit_bins1=6)

info = s_cf.plot_density(
         c2, c1, 'cross', col='mass',
         xscale='log', yscale='log', add_err=False,
         add_fit=True, fit_bins1=6, fit_log=True)

#print(info)
#info['ax'].set_xlabel('$halo_mass$')
#info['ax'].set_ylabel('LAMBSTAR')
#info['ax'].set_ylim(0,100)
#info['ax'].set_xlim(10**13.5,10**15)
#info['ax'].set_title(matching) 
plt.savefig(outpath+'scaling_log_mass_density_metrics.png', bbox_inches='tight')
plt.close()
#sys.exit()

#my mass-richness plots
#c_merged = c_merged[c_merged.data['richness']>25]
#c_merged = c_merged[c_merged.data['log_m200c']>13.5]
#M200c vs halo_mas
plt.figure()
plt.scatter(c_merged.data['log_mass'], np.log(c_merged.data['richness']), marker='.',color = 'blue', s=0.5, alpha=0.3, label='clusters')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([14.0, 15])
plt.ylim([1.0, 6])
plt.xlabel('log10(mass)')
plt.ylabel('log(richness)')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'mass_mass.png', bbox_inches='tight')

###extra plots for quick redshift studies
#M200c vs halo_mas
plt.figure()
plt.scatter(c_merged.data['z_cl'], c_merged.data['z_halo'], marker='.',color = 'blue', s=0.5, alpha=0.3, label='clusters')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0, 1.5])
plt.ylim([-0.3, 0.3])
plt.xlabel('z_cl')
plt.ylabel('z_cl-z_halo')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'redshift_redshift.png', bbox_inches='tight')

#cluster density versus m200c
sky_area_sq_deg = 440
mag_bins = np.linspace(10, 100, 90)
cdf1 = np.searchsorted(c1["mass"], mag_bins, sorter=c1["mass"].argsort())
cdf1 = len(c1)-cdf1
#cdf2 = np.searchsorted(c2["mass"], mag_bins, sorter=c2["mass"].argsort())
#cdf2 = len(c2)-cdf2
#cdf1 = np.searchsorted(c1["mass"], mag_bins, sorter=c1["mass"].argsort())
#cdf1 = len(c1)-cdf1
#cdf3 = np.searchsorted(c3["mass"], mag_bins, sorter=c3["mass"].argsort())
#cdf3 = len(c3)-cdf3
g1, = plt.semilogy(mag_bins, cdf1 / sky_area_sq_deg, color = 'red')
#g2, = plt.semilogy(mag_bins, cdf2 / sky_area_sq_deg, color = 'black')
#g3, = plt.semilogy(mag_bins, cdf3 / sky_area_sq_deg, color = 'blue')
plt.xlabel("alg. richness");
plt.ylabel("Density (cl/deg2)");
#plt.legend([g1,g2,g3],[label_1, label_2, label_3], title = 'cosmoDC2', loc='upper right')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.5', color='grey')
plt.savefig(outpath+"cluster_density.png", bbox_inches='tight')
plt.close()
