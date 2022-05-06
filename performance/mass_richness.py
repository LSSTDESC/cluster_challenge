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
from clevar.match import output_matched_catalog

matching_folder = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/after_matching/'
matching_folder = matching_folder + 'amico_cosmoDC2/'

##########select case
catalog1 = 'c1.fits'
catalog2 = 'c2.fits'
##########

outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/mass_richness/"
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
#print(c1.data)
#print(c2.data)

#create merged catalog
output_matched_catalog(matching_folder+catalog1, matching_folder+catalog2,matching_folder+'output_catalog.fits', c1, c2, matching_type='cross', overwrite=True)
c_merged = ClCatalog.read(matching_folder+'output_catalog.fits', 'merged',  cat1_z='cat1_z', cat1_mass = 'cat1_mass', cat2_z='cat2_z', cat2_mass = 'cat2_mass')
c_merged = c_merged[c_merged.data['cat2_mass']>10**13]
#print(c_merged)

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

#create 2D mass-richness histos
zbins = [0,0.5,0.75,1.0,1.25,1.5]
ybins = 10**np.linspace(13, 15, 20)
xbins = [0,10,20,30,50,70,100]
for i in range(0,5):
     cut1 = zbins[i]
     cut2 = zbins[i+1]
     filter1 = np.logical_and(c_merged.data['cat1_z'] > cut1, c_merged.data['cat1_z'] < cut2)
     c_merged_cut = c_merged[filter1]
     
     #halo_mass
     counts, xedges, yedges = np.histogram2d(c_merged_cut.data['cat1_mass'], c_merged_cut.data['cat2_mass'], bins=(xbins, ybins))
     fig, ax = plt.subplots()
     ax.pcolormesh(xbins, ybins, counts.T, cmap='jet')
     ax.set_yscale('log')
     ax.set_xscale('log')
     ax.set_xlabel('NGALS')
     ax.set_ylabel('halo_mass')
     fig.savefig(outpath+'mass_richness_halo_mass_2D_z_bin_'+str(i)+ '.png', bbox_inches='tight')
     
     file = open(outpath+'richness_z_bin_'+str(i)+ '.p', "wb")
     pickle.dump(counts, file)
     pickle.dump(xedges, file)
     pickle.dump(yedges, file)
     file.close()

#read back histograms
mass_mean = np.zeros((5,6))
mass_rms = np.zeros((5,6))
bin_x_ref = []
xedges_ref = []
yedges_ref = []

#loop over redshift bins
for redshift in range(4):
     #print('bin z ' + str(redshift))

     redshift_bin = str(round(redshift*0.1,1))+'-'+str(round(redshift*0.1+0.1,1))
     redshift_str = str(redshift)
     file = open(outpath + 'richness_z_bin_'+str(redshift)+ '.p', "rb" )
     my2D = pickle.load(file)
     xedges = pickle.load(file)
     yedges = pickle.load(file)
     xedges_ref = xedges
     yedges_ref = yedges
     h = np.array(my2D)
     #print(my2D)
     #print(xedges)
     #print(yedges)
     n_total = np.nansum(np.array(my2D),axis=None)
     #print(n_total)
     
     #project slices in 1D
     nbins_x = xedges.size-1
     nbins_y = yedges.size-1

     mass_halo_mean = np.empty([nbins_x])
     log_mass_halo_mean = np.empty([nbins_x])
     mass_halo_rms = np.empty([nbins_x])
     log_mass_halo_rms = np.empty([nbins_x])
     
     for ix in range(nbins_x):
          dm = np.empty([nbins_y])
          bin_y = np.empty([nbins_y])
          for iy in range(nbins_y):
               dm[iy] = h[ix,iy]
               bin_y[iy] = 0.5 * (yedges[iy] + yedges[iy+1])
          
          #print(dm)
          n_clusters = int(np.nansum(np.array(dm),axis=None))
          #plot
          plt.figure()
          #plt.hist(bin_y, range=bin_range, bins=nbins, weights=dm, histtype='step', color = 'black')
          plt.hist(bin_y, bins=yedges, weights=dm, histtype='step', color = 'black')
          plt.xscale('log')
          #try to fit
          #1 gauss
          #popt, pcov = curve_fit(gauss, xdata=bin_y, ydata=dm, p0=[50, 1.0, 0.05])
          #print(popt)
          # 2 gauss
          #popt, pcov = curve_fit(sum_2gauss, xdata=bin_y, ydata=dm, p0=[50, 0.5, 0.05, 50, 1.5, 0.05])
          #print(popt)
          #x = np.linspace(-0.25, 2.1, 1000)
          #fit = gauss(x, popt[0], popt[1], popt[2])
          #fit = sum_2gauss(x, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
          #gauss1 = gauss(x, popt[0], popt[1], popt[2])
          #gauss2 = gauss(x, popt[3], popt[4], popt[5])
          #plt.plot(x, gauss1, color='blue', linewidth=2.0)
          #plt.plot(x, gauss2, color='red', linewidth=2.0)
          #plt.plot(x, fit, color='red', linewidth=2.0)
          #plt.xlabel('NGALS')
          plt.xlabel('m_halo')
          #mag_bin = str(round(xedges[ix],1)) + '-' + str(round(xedges[ix+1],1))
          #plt.title('z ' + redshift_bin  + ' / mag_i ' + mag_bin)
          #plt.title('z ' + redshift_bin  + ' mag_z '+ str(xedges[bin_0]) + '-' + str(xedges[bin_0+5]))
          richness_bin = str(round(xedges[ix],1)) + '-' + str(round(xedges[ix+1],1))
          plt.title('mass_richness_bin_' + richness_bin + ': ' + str(round(n_clusters,0)) + ' clusters')
          mass_halo_mean[ix] = np.nansum(dm*bin_y) / n_clusters
          #mass_halo_mean[ix] = np.mean(dm*bin_y)
          #print(dm*bin_y)
          log_mass_halo_mean[ix] = np.log10(mass_halo_mean[ix])
          #RMS
          dx2=0
          n2=0
          for k in range(len(bin_y)):
               dx2 += dm[k]*(bin_y[k]-mass_halo_mean[ix])**2
               n2 += dm[k]
          mass_halo_rms[ix] = np.sqrt(dx2/n2)
          #mass_halo_rms[ix] = np.std(dm*bin_y)
          log_mass_halo_rms[ix] = 0.435*mass_halo_rms[ix]/mass_halo_mean[ix]
          plt.axvline(mass_halo_mean[ix],c='red')
          plt.axvline(mass_halo_mean[ix]+mass_halo_rms[ix],c='red',ls='--')
          plt.axvline(mass_halo_mean[ix]-mass_halo_rms[ix],c='red',ls='--')
          plt.savefig(outpath+'mass_in_richness_bin_' + richness_bin + '.png', bbox_inches='tight')
          plt.close()
          #print(log_mass_halo_mean[ix])
          #print(log_mass_halo_rms[ix])
          
     #summary plot
     bin_x = np.empty([nbins_x])
     for ix in range(nbins_x):
          bin_x[ix] = 0.5 * (xedges[ix] + xedges[ix+1])
     bin_x_ref = bin_x
     plt.figure()
     plt.scatter(bin_x, mass_halo_mean, label="<halo_mass>", color= "red",marker= ".", s=30)
     plt.xlim(xedges[0], xedges[nbins_x])
     #plt.ylim(yedges[0], yedges[nbins_y])
     plt.ylim(10**13, 10**15)
     #fit
     #popt, pcov = curve_fit(power, xdata=bin_x, ydata=mass_halo_mean, p0=[10**13, 0.5])
     #print(popt)
     #x = np.linspace(xedges[0], xedges[nbins_x], 1000)
     #fit = power(x, popt[0], popt[1])
     #plt.plot(x, fit, color='blue', linewidth=2.0)
     plt.xscale('log')
     plt.yscale('log')
     plt.xlabel('NGALS')
     plt.ylabel('halo_mass')
     plt.title('i-band')
     plt.legend()
     plt.savefig(outpath+"average.png", bbox_inches='tight')
     plt.close()
     #save mean
     mass_mean[redshift]=log_mass_halo_mean
     mass_rms[redshift]=log_mass_halo_rms
          
#summary of summary plots for MEAN with errors 
plt.figure()
plt.scatter(bin_x_ref, mass_mean[0], label="0-0.5", color= "black",marker= ".", s=30)
plt.plot(bin_x_ref, mass_mean[0], color= "black", linewidth=1)
#try with errors
#plt.errorbar(bin_x_ref, mass_mean[0], mass_rms[0], label="0-0.5", color= "black",marker= ".")# s=30)
#plt.plot(bin_x_ref, mass_mean[0], color= "black", linewidth=1)
plt.errorbar(bin_x_ref, mass_mean[1], mass_rms[1], label="0.5-0.75", color= "red",marker= ".")#, s=30)
plt.plot(bin_x_ref, mass_mean[1], color= "red", linewidth=1)
plt.errorbar(bin_x_ref, mass_mean[2], label="0.75-1.0", color= "blue",marker= ".")#, s=30)
plt.plot(bin_x_ref, mass_mean[2], color= "blue", linewidth=1)
plt.errorbar(bin_x_ref, mass_mean[3], label="1.0-1.25", color= "green",marker= ".")#, s=30)
plt.plot(bin_x_ref, mass_mean[3], color= "green", linewidth=1)
plt.errorbar(bin_x_ref, mass_mean[4], label="1.25-1.5", color= "purple",marker= ".")#, s=30)
plt.plot(bin_x_ref, mass_mean[4], color= "purple", linewidth=1)
plt.xlim(xedges_ref[0], xedges_ref[nbins_x])
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
#plt.ylim(yedges[0], yedges[nbins_y])
#plt.ylim(10**13, 10**15)
plt.ylim(13.8, 15)
#fit
#popt, pcov = curve_fit(power, xdata=bin_x, ydata=mass_halo_mean, p0=[10**13, 0.5])
#print(popt)
#x = np.linspace(xedges[0], xedges[nbins_x], 1000)
#fit = power(x, popt[0], popt[1])
#plt.plot(x, fit, color='blue', linewidth=2.0)
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('richness')
#plt.xlabel('NGALS_i*')
plt.ylabel('log(M200c)')
#plt.ylabel('log(halo_mass)')
plt.title('redMaPPer')
plt.legend()
plt.savefig(outpath+"average_with_errors.png", bbox_inches='tight')
plt.close()

#summary of summary plots for MEAN
#print(mass_mean)
#print(mass_rms)
plt.figure()
plt.scatter(bin_x_ref, mass_mean[0], label="0-0.5", color= "black",marker= ".", s=30)
plt.plot(bin_x_ref, mass_mean[0], color= "black", linewidth=1)
#try with errors
#plt.errorbar(bin_x_ref, mass_mean[0], mass_rms[0], label="0-0.5", color= "black",marker= ".")# s=30)
#plt.plot(bin_x_ref, mass_mean[0], color= "black", linewidth=1)
plt.scatter(bin_x_ref, mass_mean[1], label="0.5-0.75", color= "red",marker= ".", s=30)
plt.plot(bin_x_ref, mass_mean[1], color= "red", linewidth=1)
plt.scatter(bin_x_ref, mass_mean[2], label="0.75-1.0", color= "blue",marker= ".", s=30)
plt.plot(bin_x_ref, mass_mean[2], color= "blue", linewidth=1)
plt.scatter(bin_x_ref, mass_mean[3], label="1.0-1.25", color= "green",marker= ".", s=30)
plt.plot(bin_x_ref, mass_mean[3], color= "green", linewidth=1)
plt.scatter(bin_x_ref, mass_mean[4], label="1.25-1.5", color= "purple",marker= ".", s=30)
plt.plot(bin_x_ref, mass_mean[4], color= "purple", linewidth=1)
plt.xlim(xedges_ref[0], xedges_ref[nbins_x])
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
#plt.ylim(yedges[0], yedges[nbins_y])
#plt.ylim(10**13, 10**15)
plt.ylim(13.8, 15)
#fit
#popt, pcov = curve_fit(power, xdata=bin_x, ydata=mass_halo_mean, p0=[10**13, 0.5])
#print(popt)
#x = np.linspace(xedges[0], xedges[nbins_x], 1000)
#fit = power(x, popt[0], popt[1])
#plt.plot(x, fit, color='blue', linewidth=2.0)
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('richness')
#plt.xlabel('NGALS_i*')
plt.ylabel('log(M200c)')
#plt.ylabel('log(halo_mass)')
plt.title('redMaPPer')
plt.legend()
plt.savefig(outpath+"average.png", bbox_inches='tight')
plt.close() 

#summary of summary plots for RMS
plt.figure()
plt.scatter(bin_x_ref, mass_rms[0], label="0-0.5", color= "black",marker= ".", s=30)
plt.plot(bin_x_ref, mass_rms[0], color= "black", linewidth=1)
#try with errors
#plt.errorbar(bin_x_ref, mass_rms[0], mass_rms[0], label="0-0.5", color= "black",marker= ".")# s=30)
#plt.plot(bin_x_ref, mass_rms[0], color= "black", linewidth=1)
plt.scatter(bin_x_ref, mass_rms[1], label="0.5-0.75", color= "red",marker= ".", s=30)
plt.plot(bin_x_ref, mass_rms[1], color= "red", linewidth=1)
plt.scatter(bin_x_ref, mass_rms[2], label="0.75-1.0", color= "blue",marker= ".", s=30)
plt.plot(bin_x_ref, mass_rms[2], color= "blue", linewidth=1)
plt.scatter(bin_x_ref, mass_rms[3], label="1.0-1.25", color= "green",marker= ".", s=30)
plt.plot(bin_x_ref, mass_rms[3], color= "green", linewidth=1)
plt.scatter(bin_x_ref, mass_rms[4], label="1.25-1.5", color= "purple",marker= ".", s=30)
plt.plot(bin_x_ref, mass_rms[4], color= "purple", linewidth=1)
plt.xlim(xedges_ref[0], xedges_ref[nbins_x])
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
#plt.ylim(yedges[0], yedges[nbins_y])
#plt.ylim(10**13, 10**15)
plt.ylim(0, 0.5)
#fit
#popt, pcov = curve_fit(power, xdata=bin_x, ydata=mass_halo_mean, p0=[10**13, 0.5])
#print(popt)
#x = np.linspace(xedges[0], xedges[nbins_x], 1000)
#fit = power(x, popt[0], popt[1])
#plt.plot(x, fit, color='blue', linewidth=2.0)
plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('NGALS_i')
#plt.xlabel('NGALS_i*')
plt.xlabel('richness')
plt.ylabel('RMS(log(M200c))')
#plt.ylabel('log(halo_mass)')
#plt.title('i-band')
plt.legend()
plt.savefig(outpath+"average_rms.png", bbox_inches='tight')

#likelihood computation
def log_likelihood(theta, x, y, yerr):
     m0, f , g = theta
     model = m0 + f * np.log10(x[:,0]/r0) + g * np.log10((1+x[:,1])/(1+z0))
     #model = mass_richness(m0, x1, x2, f, g
     #model = m0 + f * np.log10(x/r0)
     #model = f * x + m0;
     sigma2 = yerr ** 2
     return -0.5 * np.sum((y - model) ** 2 / sigma2 )

print('Start log_likelihood minimization')
nll = lambda *args: -log_likelihood(*args)
initial = np.array([14.0, 1.0, -0.3])
#initial = np.array([14.0, 1.0])

#define x, y and yerr
##zbins = [0,0.5,0.75,1.0,1.25,1.5]
##ybins = 10**np.linspace(13, 15, 20)
##xbins = [20,25,30,40,50,60,100] 
rbins = np.array([22.5, 27.5, 35, 45., 55., 80.])
#zbins = np.array([0.35,0.65,0.85,1.1,1.35])
zbins = np.array([0.25,0.625,0.875,1.125])
#r0 = np.linspace(10, 100, 1000)
#z0 = np.linspace(0, 1.5, 1000)
n_points = len(rbins)*len(zbins)
x0 = np.empty((n_points,2))
#x = np.empty(n_points)
x0_2 = np.empty(n_points)
y = np.empty(n_points)
y_err = np.empty(n_points)
print('number of points = ' + str(n_points))
for i in range(len(rbins)-1):
     for j in range(len(zbins)):
          print(str(i) + ' ' + str(j))
          k = i * (len(rbins)-1) + j
          print(str(k))
          x0[k] =  np.array([rbins[i],zbins[j]])
          #x[k] = rbins[i]
          x0_2[k] = zbins[j]
          y[k] = mass_mean[j][i]
          #y_err[k] = 0.001*y[k] #mass_rms[j][i]
          y_err[k] = mass_rms[j][i]
          if y_err[k]==0:
               y_err[k] = 0.2*y[k]
print('-------')
print(x0)
#print(x0_2)
print(y)
print(y_err)
#minimization
soln = minimize(nll, initial, args=(x0, y, y_err))
m_ml, f_ml, g_ml= soln.x
print('likelihood results')
print(m_ml)
print(f_ml)
print(g_ml)

print('*********emcee**********')
initial_guess = np.array([14, 1, 0.02])
nwalkers = 100
p0 = np.random.randn(nwalkers, initial_guess.size)*0.1+initial_guess
sampler = emcee.EnsembleSampler(
                 nwalkers, initial_guess.size,
                 log_likelihood, args=[x0, y, y_err])
sampler.run_mcmc(p0, 100)
#chains = sampler.get_chain()
chains = [sampler.get_chain()[:,:,i].flatten() for i in range(initial_guess.size)]
#check chains
plt.figure()
plt.plot(chains[0])
plt.savefig(outpath+"chain0.png")
plt.close()
plt.figure()
plt.plot(chains[1])
plt.savefig(outpath+"chain1.png")
plt.close()
#threshold at nchains * nwalkers?
mean0 = np.mean(chains[0][10000:])
mean1 = np.mean(chains[1][10000:])
mean2 = np.mean(chains[2][10000:])

print(str(mean0) + '+/-' + str(mean_err0))
print(str(mean1) + '+/-' + str(mean_err1))
print(str(mean2) + '+/-' + str(mean_err2))
flat_chains = sampler.get_chain(discard=100, flat=True)
labels = ["m0", "F", "G"]
fig = corner.corner(
             flat_chains,
             labels=labels,
             color = 'blue',
             smooth = True,
             plot_datapoints = False,
             plot_density = False,
             plot_contours = True,
             fill_contours = True,
             levels = (0.68 , 0.95),
             quantiles=(0.16,0.84),
             show_titles=True
             );
#truths=[14.0, 1.0, 0.2]
fig
fig.savefig(outpath+"posteriors_2D_fit.png")

sys.exit()
