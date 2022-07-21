#!/usr/bin/env python
# coding: utf-8
# Author: T. Guillemin
#######################################################
# Adapted from Constantin's Payerne notebook: Unbinned_Mass_richness_relation.ipynb
#######################################################

###import
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
import sys
import os
import shutil
import pickle
import h5py
import pandas as pd
from scipy import optimize
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.integrate import simps
from scipy.stats import norm
from scipy.stats import lognorm
from scipy import optimize
from scipy.optimize import minimize
from numpy import exp, linspace, random
import emcee
import corner

###plot style
plt.rcParams['figure.figsize'] = [9.5, 6]
plt.rcParams.update({'font.size': 18})

###inputs        
outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/selection_function/global_fit/"
if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)

#inpath = "/sps/lsst/users/ccombet/SkySim5000/hdf5/"
inpath = "/sps/lsst/users/tguillem/catalogs/SkySim5000/hdf5/"

#read hdf5 files
with pd.HDFStore(os.path.join(inpath,f'skysim_halos_z=0-1.20_mfof_gt_1.00e+13_small.hdf5')) as store:
     halo_data = store['skysim']
     halo_metadata = store.get_storer('skysim').attrs.metadata
#print(halo_data['baseDC2/sod_halo_mass'])

#rename
halo_data.rename(columns={'baseDC2/sod_halo_mass': 'M200c', 'richness': 'NGALS', 'richness_i': 'NGALS_i', 'richness_z': 'NGALS_z'}, inplace=True)
halo_data = halo_data[halo_data['M200c']>10*14.5]
#print(halo_data)
#fix M200c
halo_data['M200c'] = halo_data['M200c']/0.71

#get mean values
z_0 = np.mean(halo_data['redshift'])
logM_0 = np.mean(np.log10(halo_data['M200c']))

###fit functions
def gaussian(x, mu, sigma):
     return 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

def log_normal(x, mu, sigma):
     return 1/(x*np.sqrt(2*np.pi*sigma**2))*np.exp(-(np.log(x)-mu)**2/(2*sigma**2))

def mu_lnlambda(redshift, logM, mu0, G_z_mu, G_logM_mu):
         return mu0 + G_z_mu*(1+redshift)/(1+z_0) + G_logM_mu*(logM)/logM_0
    
def sigma_lnlambda(redshift, logM, sigma0, F_z_sigma, F_logM_sigma):
         return sigma0 + F_z_sigma*(1+redshift)/(1+z_0) + F_logM_sigma*(logM)/logM_0

def P(data, theta):
     r"""probability of ln(lambda) knowing the value of theta"""
     mu0, G_z_mu, G_logM_mu, sigma0, F_z_sigma, F_logM_sigma = theta
     lnlambda, redshift, logM = data
     mu = mu_lnlambda(redshift, logM, mu0, G_z_mu, G_logM_mu)
     sigma = sigma_lnlambda(redshift, logM, sigma0, F_z_sigma, F_logM_sigma)
     return gaussian(lnlambda, mu, sigma)

#input preparation for the fit
logM = np.log10(halo_data['M200c'])
lnlambda= np.log(halo_data['NGALS'])
redshift = halo_data['redshift']
data_mcmc = [lnlambda, redshift, logM]

#likelihood definition
def lnL(theta):
     p = P(data_mcmc, theta)
     return np.sum(np.log(p))

#MCMC
n_chains = 10000
initial = [0.5,25.0,-25.0,-0.10,-2.5,2.5]
pos = initial + 0.001 * np.random.randn(70, len(initial))
nwalkers, ndim = pos.shape
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnL)
sampler.run_mcmc(pos, n_chains, progress=False); 

fig, axes = plt.subplots(6, figsize=(10, 7), sharex=True)
labels = [r'$\mu_0})$', r'$G_z^\mu$', r'$G_logM^\mu$', r'$\sigma_0$', r'$F_z^\sigma$', r'$F_logM^\sigma$']
samples = sampler.get_chain()
for i in range(6):
         ax = axes[i]
         ax.plot(samples[:, :, i], "k", alpha=0.3)
         ax.set_xlim(0, len(samples))
         ax.set_ylabel(labels[i])
         ax.yaxis.set_label_coords(-0.1, 0.5)
         
flat_samples = sampler.get_chain(discard=10, thin=1,flat=True)

def mu_sigma_parameters(samples = None, labels = None):
     mu = {labels[i] : np.mean(samples[:,i]) for i in range(len(labels))}
     sigma = {labels[i] : np.std(samples[:,i]) for i in range(len(labels))}
     return mu, sigma

sigma1 = 1. - np.exp(-(1./1.)**2/2.)
sigma2 = 1. - np.exp(-(2./1.)**2/2.)
sigma3 = 1. - np.exp(-(3./1.)**2/2.)

bins = 20
fig, axs = plt.subplots(len(labels), len(labels), figsize = (10,10))
for i in range(len(labels)):
     for j in range(len(labels)):
          axs[i,j].tick_params(axis='both', which = 'major', labelsize= 9)
          fig = corner.corner(
               flat_samples,
               bins=bins, levels=(sigma1, sigma2, sigma3),
               fig = fig,
               smooth1d=False,smooth=False,plot_datapoints=True,
               fill_contours=True,labels = labels,
               color='deepskyblue',
               label_kwargs={"fontsize": 20},
               use_math_text=True,
               plot_density=False,
               max_n_ticks = 5,
          ); 
plt.savefig(outpath+'posteriors', bbox_inches='tight')
sys.exit()

mu_p, sigma_p = mu_sigma_parameters(samples = flat_samples, labels = labels)

richness = np.linspace(20, 200, 100)
c = ['y', 'orange', 'r', 'b',]
fig, ax = plt.subplots(1, 2, figsize = (14,5))
ax[0].scatter(logrichness, logm200c, alpha = .2)
for i, z in enumerate([.1, .2, .5, 1]):
     mu = mu_logM_lambda(z, np.log10(richness), mu_p[r'$\log_{10}(M_{\rm 200c,0})$'], mu_p[r'$G_z^\mu$'], mu_p[r'$G_\lambda^\mu$'])
     sigma = sigma_logM_lambda(z, np.log10(richness), mu_p[r'$\sigma_{\log_{10}M}$'], mu_p[r'$F_z^\sigma$'], mu_p[r'$F_\lambda^\sigma$'])
     ax[0].plot(np.log10(richness), mu, c = c[i] , label = f'z = {z:.2f}', linewidth = 3)
     ax[0].legend()
     ax[1].plot(np.log10(richness), sigma, c = c[i] , label = f'z = {z:.2f}', linewidth = 3)
     ax[1].legend()
     ax[0].set_ylabel(r'$\log_{10}M$', fontsize = 20)
     ax[1].set_ylabel(r'$\sigma (\log_{10}M)$', fontsize = 20)
for i in range(2):
     ax[i].set_xlabel(r'$\log_{10}\lambda$', fontsize = 20)
                                                     





sys.exit()









#richness plots in (m,z) bins
zbins = [0,0.5,0.75,1.0,1.2]
ybins = 10**np.linspace(13, 14.4, 8)
n_z = len(zbins)-1
n_y = len(ybins)-1
#print(ybins)
#for j in range(0,n_y):
     #print(ybins[j])
a_mu = np.zeros((n_z,n_y))
a_sigma = np.zeros((n_z,n_y))
a_z = np.zeros((n_z,n_y))
a_mass = np.zeros((n_z,n_y))
for i in range(0,n_z):
     cut1 = zbins[i]
     cut2 = zbins[i+1]
     filter1 = np.logical_and(halo_data['redshift'] > cut1, halo_data['redshift'] < cut2)
     halos_1 = halo_data[filter1]
     for j in range(0,n_y):
          cut3 = ybins[j]
          cut4 = ybins[j+1]
          #print(cut3)
          filter2 = np.logical_and(halos_1['M200c'] > cut3, halos_1['M200c'] < cut4)
          halos = halos_1[filter2]
          #print(halos)
          #get means of z distribution and of M200c definition
          a_z[i][j] = np.mean(halos['redshift'])
          a_mass[i][j] = np.mean(np.log10(halos['M200c']))
          #richness
          nbins = 40
          bin_range = [0,200]
          plt.figure()
          plt.hist(halos['NGALS_i'], range=bin_range, bins=nbins, label='Halos', histtype='step', color = 'black', density=True)#, stacked=True)
          #plt.hist(halos['NGALS'], range=bin_range, bins=nbins, label='NGALS', histtype='step', color = 'black')
          #plt.hist(halos['NGALS_i'], range=bin_range, bins=nbins, label='NGALS_i', histtype='step', color = 'red')
          #plt.hist(halos['NGALS_z'], range=bin_range, bins=nbins, label='NGALS_z', histtype='step', color = 'blue')
          plt.xlabel("r");
          plt.ylabel("P(r|m,z)")
          #plt.xscale('log')
          #plt.yscale('log')
          plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
          plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
          f_cut1=round(cut1,1)
          f_cut2=round(cut2,1)
          f_cut3=round(np.log10(cut3),1)
          f_cut4=round(np.log10(cut4),1)
          plt.title('cosmoDC2: m200c '+str(f_cut3)+'-'+str(f_cut4) + ' / z '+str(f_cut1)+'-'+str(f_cut2))
          #plt.legend(title = '', loc='upper right')
          #test fit
          print('***********FIT**************')
          xbins = np.linspace(bin_range[0], bin_range[1], nbins+1)
          counts, xedges = np.histogram(halos['NGALS'],density=True,bins=xbins)
          #print(counts)
          #print(xedges)
          bin_x = np.empty([nbins])
          for i_bin in range(nbins):
               bin_x[i_bin] = 0.5 * (xedges[i_bin] + xedges[i_bin+1])
          #print(len(bin_x))
          #print(len(counts))
          #dr = np.empty([nbins])
          #bin_x = np.empty([nbins])
          ###for i_bin in range(nbins):
          ###     dr[i_bin] = h[ix,iy]
          ###bin_y[iy] = 0.5 * (yedges[iy] + yedges[iy+1])
          #popt, pcov = curve_fit(gauss, xdata=bin_x, ydata=counts, p0=[0.05, 20, 5])
          popt, pcov = curve_fit(log_normal, xdata=bin_x, ydata=counts, p0=[8, 3])
          print(popt)
          f_popt=np.around(popt,decimals=2)
          x = np.linspace(0.1, 200, 2000) 
          #gauss1 = gauss(x, popt[0], popt[1], popt[2])
          #plt.plot(x, gauss1, color='blue', linewidth=2.0)
          log_normal_1 = log_normal(x, popt[0], popt[1])
          plt.plot(x, log_normal_1, color='blue', linewidth=2.0,label="Fit "+str(f_popt))
          plt.legend()
          #mu, std = lognorm.fit(halos['NGALS']) 
          #print(mu)
          #print(std)
          #xmin, xmax = plt.xlim()

          #x = np.linspace(xmin, xmax, 100)
          #p = norm.pdf(x, mu, std)
          #plt.plot(x, p, 'k', linewidth=2)
          plt.savefig(outpath+'richness_redshift_bin_'+str(i)+'_mass_bin_'+str(j)+'.png')
          plt.close()
          #save fit parameters
          a_mu[i][j]=popt[0]
          a_sigma[i][j]=abs(popt[1])

#parametrize mu and sigma
#print(a_mu)
#print(a_sigma)
print(a_z)
print(a_mass)
sys.exit()

#summary plots for a_mu versus log(m)
bin_x1=np.empty([n_z]) 
for ix1 in range(n_z):
     bin_x1[ix1] = 0.5 * (zbins[ix1] + zbins[ix1+1])
print(bin_x1)
bin_x2=np.empty([n_y]) 
for ix2 in range(n_y):
     bin_x2[ix2] = np.log10(0.5 * (ybins[ix2] + ybins[ix2+1]))
print(bin_x2)
#summary plot for a_mu versus log(m)
plt.figure()
plt.scatter(bin_x2, a_mu[0,:], label="0-0.5", color= "black",marker= ".", s=30)
plt.plot(bin_x2, a_mu[0,:], color= "black", linewidth=1)
plt.scatter(bin_x2, a_mu[1,:], label="0.5-0.75", color= "red",marker= ".", s=30)
plt.plot(bin_x2, a_mu[1,:], color= "red", linewidth=1)
plt.scatter(bin_x2, a_mu[2,:], label="0.75-1.0", color= "blue",marker= ".", s=30)
plt.plot(bin_x2, a_mu[2,:], color= "blue", linewidth=1)
plt.scatter(bin_x2, a_mu[3,:], label="1.0-1.2", color= "green",marker= ".", s=30)
plt.plot(bin_x2, a_mu[3,:], color= "green", linewidth=1)
plt.xlabel('log(m200c)')
plt.legend()
plt.savefig(outpath+"mean_vs_mass.png", bbox_inches='tight')
plt.close()
#summary plot for a_sigma versus log(m)
plt.figure()
plt.scatter(bin_x2, a_sigma[0,:], label="0-0.5", color= "black",marker= ".", s=30)
plt.plot(bin_x2, a_sigma[0,:], color= "black", linewidth=1)
plt.scatter(bin_x2, a_sigma[1,:], label="0.5-0.75", color= "red",marker= ".", s=30)
plt.plot(bin_x2, a_sigma[1,:], color= "red", linewidth=1)
plt.scatter(bin_x2, a_sigma[2,:], label="0.75-1.0", color= "blue",marker= ".", s=30)
plt.plot(bin_x2, a_sigma[2,:], color= "blue", linewidth=1)
plt.scatter(bin_x2, a_sigma[3,:], label="1.0-1.2", color= "green",marker= ".", s=30)
plt.plot(bin_x2, a_sigma[3,:], color= "green", linewidth=1)
plt.xlabel('m200c')
plt.legend()
plt.savefig(outpath+"sigma_vs_mass.png", bbox_inches='tight')
plt.close()
#summary plot for a_mu versus z
plt.figure()
plt.scatter(bin_x1, a_mu[:,0], label="mass bin 1", color= "black",marker= ".", s=30)
plt.plot(bin_x1, a_mu[:,0], color= "black", linewidth=1)
plt.scatter(bin_x1, a_mu[:,1], label="mass bin 2", color= "red",marker= ".", s=30)
plt.plot(bin_x1, a_mu[:,1], color= "red", linewidth=1)
plt.scatter(bin_x1, a_mu[:,2], label="mass bin 3", color= "blue",marker= ".", s=30)
plt.plot(bin_x1, a_mu[:,2], color= "blue", linewidth=1)
plt.scatter(bin_x1, a_mu[:,3], label="mass bin 4", color= "green",marker= ".", s=30)
plt.plot(bin_x1, a_mu[:,3], color= "green", linewidth=1)
plt.scatter(bin_x1, a_mu[:,4], label="mass bin 5", color= "purple",marker= ".", s=30)
plt.plot(bin_x1, a_mu[:,4], color= "purple", linewidth=1)
#plt.scatter(bin_x1, a_mu[:,5], label="mass bin 6", color= "orange",marker= ".", s=30)
#plt.plot(bin_x1, a_mu[:,5], color= "orange", linewidth=1)
#plt.scatter(bin_x1, a_mu[:,6], label="mass bin 7", color= "brown",marker= ".", s=30)
#plt.plot(bin_x1, a_mu[:,6], color= "brown", linewidth=1)
plt.xlabel('z')
#plt.legend()
plt.savefig(outpath+"mean_vs_z.png", bbox_inches='tight')
plt.close()
#summary plot for a_sigma versus z
plt.figure()
plt.scatter(bin_x1, a_sigma[:,0], label="mass bin 1", color= "black",marker= ".", s=30)
plt.plot(bin_x1, a_sigma[:,0], color= "black", linewidth=1)
plt.scatter(bin_x1, a_sigma[:,1], label="mass bin 2", color= "red",marker= ".", s=30)
plt.plot(bin_x1, a_sigma[:,1], color= "red", linewidth=1)
plt.scatter(bin_x1, a_sigma[:,2], label="mass bin 3", color= "blue",marker= ".", s=30)
plt.plot(bin_x1, a_sigma[:,2], color= "blue", linewidth=1)
plt.scatter(bin_x1, a_sigma[:,3], label="mass bin 4", color= "green",marker= ".", s=30)
plt.plot(bin_x1, a_sigma[:,3], color= "green", linewidth=1)
plt.scatter(bin_x1, a_sigma[:,4], label="mass bin 5", color= "purple",marker= ".", s=30)
plt.plot(bin_x1, a_sigma[:,4], color= "purple", linewidth=1)
#plt.scatter(bin_x1, a_sigma[:,5], label="mass bin 6", color= "orange",marker= ".", s=30)
#plt.plot(bin_x1, a_sigma[:,5], color= "orange", linewidth=1)
#plt.scatter(bin_x1, a_sigma[:,6], label="mass bin 7", color= "brown",marker= ".", s=30)
#plt.plot(bin_x1, a_sigma[:,6], color= "brown", linewidth=1)
plt.xlabel('z')
#plt.legend()
plt.savefig(outpath+"sigma_vs_z.png", bbox_inches='tight')
plt.close()

#likelihood computation
def log_likelihood(theta, x, y, yerr):
          a, b, c = theta
          #model =  a*np.log10((1+x[:,0])/(1+0.8)) + b*x[:,1]/np.log10(13.6) + np.log10(c)
          model =  a*(1+x[:,0])/(1+0.8) + b*x[:,1]/13.5 + c
          sigma2 = yerr ** 2
          return -0.5 * np.sum((y - model) ** 2 / sigma2 )

print('Start log_likelihood minimization')
nll = lambda *args: -log_likelihood(*args)
initial = np.array([10,10,10])

#define x, y and yerr
n_points = n_z*n_y
x0 = np.empty((n_points,2))
y = np.empty(n_points)
y_s = np.empty(n_points)
y_err = np.empty(n_points)
y_s_err = np.empty(n_points)
print('number of points = ' + str(n_points))
for i in range(n_z):
     for j in range(n_y):
          print(str(i) + ' ' + str(j))
          k = i * n_y  + j
          print(str(k))
          x0[k] =  np.array([bin_x1[i],bin_x2[j]])
          y[k] = a_mu[i][j]
          y_err[k] = 0.05*y[k]
          if y_err[k]==0:
               y_err[k] = 0.2*y[k]
          y_s[k] = a_sigma[i][j]
          y_s_err[k] = 0.05*y_s[k]
          if y_s_err[k]==0:
               y_s_err[k] = 0.2*y_s[k]     
print('-------')
print(x0)
print(y)
print(y_err)
print(y)
print(y_err)

#likelihood minimization minimization
#mean
soln = minimize(nll, initial, args=(x0, y, y_err))
a_ml, b_ml, c_ml= soln.x
print('likelihood results')
print(a_ml)
print(b_ml)
print(c_ml) 

def f_a_mu(x,y,a,b,c):
     #return a*np.log10((1+x)/(1+0.8)) + b*y/np.log(13.6) + np.log(c) 
     return a*(1+x)/(1+0.8) + b*y/13.5 + c
print('0.65,13.9')
print(f_a_mu(0.65,13.9,a_ml,b_ml,c_ml))
print('0.90,13.2')
print(f_a_mu(0.90,13.2,a_ml,b_ml,c_ml))

#sigma
soln_s = minimize(nll, initial, args=(x0, y_s, y_s_err))
a_s_ml, b_s_ml, c_s_ml= soln_s.x
print('likelihood results')
print(a_s_ml)
print(b_s_ml)
print(c_s_ml) 

print('Closure tests')
def f_a_sigma(x,y,a,b,c):
     return a*(1+x)/(1+0.8) + b*y/13.5 + c
print('0.65,13.9')
print(f_a_sigma(0.65,13.9,a_s_ml,b_s_ml,c_s_ml))
print('0.90,13.2')
print(f_a_sigma(0.90,13.2,a_s_ml,b_s_ml,c_s_ml))

#comparison of fit and parametrization
for i in range(0,n_z):
     cut1 = zbins[i]
     cut2 = zbins[i+1] 
     for j in range(0,n_y):
          cut3 = ybins[j]
          cut4 = ybins[j+1] 
          x = np.linspace(0.1, 200, 2000) 
          log_normal_1 = log_normal(x, a_mu[i][j], a_sigma[i][j])
          plt.plot(x, log_normal_1, color='blue', linewidth=2.0,label="Fit")
          k = i * n_y  + j
          mu = f_a_mu(x0[k][0],x0[k][1],a_ml,b_ml,c_ml)
          sigma = f_a_sigma(x0[k][0],x0[k][1],a_s_ml,b_s_ml,c_s_ml)
          log_normal_2 = log_normal(x, mu, sigma)
          plt.plot(x, log_normal_2, color='red', linewidth=2.0,label="Param")
          plt.legend()
          f_cut1=round(cut1,1)
          f_cut2=round(cut2,1)
          f_cut3=round(np.log10(cut3),1)
          f_cut4=round(np.log10(cut4),1) 
          plt.title('cosmoDC2: m200c '+str(f_cut3)+'-'+str(f_cut4) + ' / z '+str(f_cut1)+'-'+str(f_cut2))
          plt.savefig(outpath+'richness_redshift_bin_'+str(i)+'_mass_bin_'+str(j)+'_closure.png')
          plt.close()

print('*********emcee**********')
initial_guess = np.array([0.5, 25, -24])
nwalkers = 100
p0 = np.random.randn(nwalkers, initial_guess.size)*0.1+initial_guess
sampler = emcee.EnsembleSampler(
                      nwalkers, initial_guess.size,
                      log_likelihood, args=[x0, y, y_err])
sampler.run_mcmc(p0, 10000)
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
plt.figure()
plt.plot(chains[2])
plt.savefig(outpath+"chain2.png")
plt.close()
#threshold at nchains * nwalkers?
mean0 = np.mean(chains[0][10000:])
mean1 = np.mean(chains[1][10000:])
mean2 = np.mean(chains[2][10000:]) 

flat_chains = sampler.get_chain(discard=100, flat=True)
labels = [r'$a_{\mu}$', r'$b_{\mu}$', r'$c_{\mu}$']
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
                  show_titles=False
                  );
fig.savefig(outpath+"posteriors_2D_fit_mean.png")

#sigma
initial_guess = np.array([0.0, -2.0, 2.0])
nwalkers = 100
p0 = np.random.randn(nwalkers, initial_guess.size)*0.1+initial_guess
sampler = emcee.EnsembleSampler(
                      nwalkers, initial_guess.size,
                      log_likelihood, args=[x0, y_s, y_s_err])
sampler.run_mcmc(p0, 10000)
#chains = sampler.get_chain()
chains = [sampler.get_chain()[:,:,i].flatten() for i in range(initial_guess.size)]
#check chains
plt.figure()
plt.plot(chains[0])
plt.savefig(outpath+"chain0_s.png")
plt.close()
plt.figure()
plt.plot(chains[1])
plt.savefig(outpath+"chain1_s.png")
plt.close()
plt.figure()
plt.plot(chains[2])
plt.savefig(outpath+"chain2_s.png")
plt.close()
#threshold at nchains * nwalkers?
mean0 = np.mean(chains[0][10000:])
mean1 = np.mean(chains[1][10000:])
mean2 = np.mean(chains[2][10000:]) 

flat_chains = sampler.get_chain(discard=100, flat=True)
labels = [r'$a_{\sigma}$', r'$b_{\sigma}$', r'$c_{\sigma}$']
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
                  show_titles=False
                  );
fig.savefig(outpath+"posteriors_2D_fit_sigma.png")

sys.exit()
