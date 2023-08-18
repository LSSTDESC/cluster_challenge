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
from clevar.match import MembershipMatch 
from clevar.match import get_matched_pairs
from clevar.match_metrics import scaling
from clevar.match_metrics import recovery
from clevar.match_metrics import distances
from clevar.match_metrics.recovery import ClCatalogFuncs as r_cf
from clevar.match import output_matched_catalog

matching_folder = '/sps/lsst/groups/clusters/redmapper_validation_project/cosmoDC2_v1.1.4/extragal/after_matching/v1/'
#matching_folder = '/sps/lsst/groups/clusters/wazp_validation_project/cosmoDC2_v1.1.4/extragal/after_matching/v0/'
#matching_folder = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/after_matching/member_snr_gt4/'

#outpath_base = '/pbs/home/t/tguillem/web/clusters/cluster_challenge/selection_function/completeness_fit/member_snr_gt4/'
outpath_base = '/pbs/home/t/tguillem/web/clusters/cluster_challenge/debug/'

#select the catalogs to match
wazp_cosmoDC2 = False
redmapper_cosmoDC2 = True
amico_cosmoDC2 = False
matching_selected = 'cross'

if wazp_cosmoDC2 == True:
     matching_folder = matching_folder + 'wazp_cosmoDC2/'
     matching = 'WaZP cosmoDC2: NGALS>0, $m_{200c}>10^{13}$'
     outpath = outpath_base + 'wazp_cosmoDC2/'
elif redmapper_cosmoDC2 == True:
     matching_folder = matching_folder + 'redmapper_cosmoDC2/'
     matching = 'redMaPPer cosmoDC2: $\Lambda>0$, $m_{200c}>10^{13}$'
     outpath = outpath_base + 'redmapper_cosmoDC2/'
elif amico_cosmoDC2 == True:
     matching_folder = matching_folder + 'amico_cosmoDC2/'
     matching = 'AMICO cosmoDC2: $m_{halo}>10^{13}$'
     outpath = outpath_base + 'amico_cosmoDC2/'
else:
     print('Catalog selection is wrong.')
     sys.exit()

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)

##########select case
catalog1 = 'c1.fits'
catalog2 = 'c2.fits'
##########

#load c1 and c2
c1 = ClCatalog.read_full(matching_folder + catalog1)
c2 = ClCatalog.read_full(matching_folder + catalog2)
#mt = MembershipMatch()
#mt.load_shared_members(c1, c2, fileprefix=matching_folder+'mem_share')
#print(c1.mt_input)
#print(c1.data)
#print(c2.data)
#sys.exit()
#restrict to z<1.15
#c1=c1[c1.data['z']<1.15]

#restrict to matched pairs
#mt1, mt2 = get_matched_pairs(c1, c2, 'cross', None, None) 

#plot style
figx=10
figy=7

#recovery_plot
#zbins = np.linspace(0.2,1.5,14)
zbins = np.linspace(0.2,1.1,12)
mbins = [10**13,10**13.5,10**14,10**14.5,10**15]
#zbins = np.linspace(0, 1.6, 9)
#mbins = np.logspace(14, 15, 5)
fig = plt.figure()#figsize=(figx,figy))
info = r_cf.plot(c2, col1='z', col2='mass', bins1=zbins, bins2=mbins, matching_type=matching_selected, legend_format=lambda x: f'10^{{{np.log10(x)}}}', lines_kwargs_list = [{'color':'black'}, {'color':'red'}, {'color':'blue'}, {'color':'purple'}])
info['ax'].set_xlabel('$z_{halo}$')
info['ax'].set_ylabel('Completeness') 
info['ax'].set_ylim(0,1.2)
info['ax'].set_xlim(0.2,1.6) 
info['ax'].set_title(matching)
plt.savefig(outpath+'recovery_plot.png', bbox_inches='tight')
plt.close(fig)

#recovery_plot_panel
mbins = np.logspace(13, 15, 10)
zbins = np.linspace(0.2,1.2,11)
plt.figure()
info = r_cf.plot_panel(c2, col1='z', col2='mass', bins1=zbins, bins2=mbins,
                       matching_type='multi_self', label_format=lambda x: f'10^{{{np.log10(x)}}}')
plt.savefig(outpath+'recovery_plot_panel.png', bbox_inches='tight')
plt.close()

#recovery_plot2D linear
#mbins = np.linspace(13, 15, 17)
#plt.figure()
#info = r_cf.plot2D(c2, col1='z', col2='log_mass', bins1=zbins, bins2=mbins,
#                   matching_type=matching_selected, plt_kwargs={'cmap':'jet'})#{'cmap':'nipy_spectra'} #jet #gnuplot
#plt.savefig(outpath+'recovery_plot2D.png', bbox_inches='tight') 

#purity
zbins = np.linspace(0.2,1.1,10)
mbins = [0,10,20,30,100]
#zbins = np.linspace(0, 1.6, 9)
#mbins = np.logspace(14, 15, 5)
fig = plt.figure()
info = r_cf.plot(c1, col1='z', col2='mass', bins1=zbins, bins2=mbins, matching_type=matching_selected, lines_kwargs_list = [{'color':'black'}, {'color':'red'}, {'color':'blue'}, {'color':'purple'}])
info['ax'].set_xlabel('$z_{cl}$')
info['ax'].set_ylabel('Purity') 
info['ax'].set_ylim(0.5,1.2)
info['ax'].set_xlim(0.2,1.6) 
info['ax'].set_title(matching)
plt.savefig(outpath+'purity_plot.png', bbox_inches='tight')
plt.close(fig)

#recovery_plot2D
plt.figure()
info = r_cf.plot2D(c1, col1='z', col2='mass', bins1=zbins, bins2=mbins,
                   matching_type=matching_selected, plt_kwargs={'cmap':'jet'})
plt.savefig(outpath+'purity_plot2D.png', bbox_inches='tight') 
#sys.exit()

#purity versus mass
zbins = np.linspace(0, 1.5, 2)
mbins = [5,10,15,20,30,50,100]
plt.figure()
info = recovery.plot(c1, matching_selected, zbins, mbins,
                     shape='line', transpose=True)
plt.savefig(outpath+'purity_richness.png', bbox_inches='tight')

###EXTRA CHECKS
#ra/dec map
plt.figure()
plt.plot(c1.data['ra'],c1.data['dec'],'rx', color = 'red', label = 'clusters')
plt.plot(c2.data['ra'],c2.data['dec'],'rx', color = 'blue', label = 'halos')
plt.xlim([60, 72])
plt.ylim([-50, -30])
plt.xlabel("ra")
plt.ylabel("dec");
plt.title(matching)
plt.legend(bbox_to_anchor = (1, 1), loc = 'upper right', prop = {'size': 15})
plt.savefig(outpath+"clusters.png", bbox_inches='tight')
plt.close()

#mass
#path_c3 = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/batch/after_matching/full/m200c/wazp_cosmoDC2/c1.fits'
#path_c4 = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/batch/after_matching/full/m200c/wazp_cosmoDC2/c2.fits'
#c3  = ClCatalog.read_full(path_c3)
#c4  = ClCatalog.read_full(path_c4)
#c4
#c1=c1[c1.data['mass']>25]
#c1=c1[c1.data['m200c']>10**14]
#c3=c3[c3.data['mass']>20]
#c3=c1[c3.data['m200c']>10**14]

#create a merged catalog for the cross-matched pairs
output_matched_catalog(matching_folder+catalog1, matching_folder+catalog2,matching_folder+'output_catalog_12.fits', c1, c2, matching_type='cross', overwrite=True)
c_merged_12 = ClCatalog.read(matching_folder+'output_catalog_12.fits', 'merged',  z_cl='cat1_z', richness = 'cat1_mass', z_halo='cat2_z', mass = 'cat2_mass', log_mass = 'cat2_log_mass', m200c = 'cat2_m200c', log_m200c = 'cat2_log_m200c',)
c_merged_12 = c_merged_12[c_merged_12.data['z_cl']<1.15]
#c_merged_12 = c_merged_12[c_merged_12.data['richness']>20]
#output_matched_catalog(path_c3, path_c4,matching_folder+'output_catalog_34.fits', c3, c4, matching_type='cross', overwrite=True)
#c_merged_34 = ClCatalog.read(matching_folder+'output_catalog_34.fits', 'merged',  z_cl='cat1_z', richness = 'cat1_mass', z_halo='cat2_z', mass = 'cat2_mass', log_mass = 'cat2_log_mass', m200c = 'cat2_m200c', log_m200c = 'cat2_m200c',)
#c_merged_34 = c_merged_34[c_merged_34.data['m200c']>10**14]
#c_merged_34 = c_merged_34[c_merged_34.data['richness']>25]
bin_range = [0,60]
nbins = 60
plt.figure()
plt.hist(c_merged_12.data['richness'], range=bin_range, bins=nbins, label='RM', histtype='step', color = 'red')
#plt.hist(c_merged_34.data['richness'], range=bin_range, bins=nbins, label='WaZP', histtype='step', color = 'black')
plt.xlabel("alg. richness");
plt.ylabel("clusters")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.legend(title = '', loc='upper right')
plt.title('cosmoDC2')
plt.savefig(outpath+'mass.png', bbox_inches='tight')
plt.close()

#redshift
bin_range = [0,1.6]
nbins = 32
plt.figure()
plt.hist(c_merged_12.data['z_cl'], range=bin_range, bins=nbins, label='RM', histtype='step', color = 'red')
#plt.hist(c_merged_34.data['z_cl'], range=bin_range, bins=nbins, label='WaZP', histtype='step', color = 'black')
#plt.hist(c3['z'], range=bin_range, bins=nbins, label=label_3, histtype='step', color = 'blue')
plt.xlabel("redshift");
plt.ylabel("clusters / 0.05 dz")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.legend(title = '', loc='upper right')
plt.title('cosmoDC2')
plt.savefig(outpath+'redshift.png', bbox_inches='tight')
plt.close()

###selection functions
#in redshift bins, check versus mass
###completeness
cluster_unmatched = c1[c1.data['mt_cross']==None]
halos_unmatched = c2[c2.data['mt_cross']==None]

#versus z
print('++++++++++++++++Completeness in redshift bins')
bin_range = [0.2,1.1]
nbins_x = 9
x_bins = np.linspace(0.2,1.1,nbins_x+1)
mbins = [10**13,10**13.5,10**14,10**14.5,10**15]
compl_z_raw = np.empty([4,nbins_x])
compl_z = np.empty([4,nbins_x])

for i in range(0,4):
     print('-----'+str(i))
     cut1 = mbins[i]
     cut2 = mbins[i+1]
     filter1 = np.logical_and(c_merged_12.data['m200c'] > cut1, c_merged_12.data['m200c'] < cut2)
     c_halos_matched = c_merged_12[filter1]
     #print(c_halos_matched)
     filter2 = np.logical_and(c2.data['m200c'] > cut1, c2.data['m200c'] < cut2)
     c_halos = c2.data[filter2]
     #print(c_halos)
     h_r_halos_matched = np.histogram(c_halos_matched['z_halo'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
     h_r_halos  = np.histogram(c_halos['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
     print(h_r_halos_matched)
     print(h_r_halos)
     compl_z_raw[i] = np.divide(h_r_halos_matched[0],h_r_halos[0],where=(h_r_halos[0]!=0))
     print(compl_z_raw[i])

#versus mass
print('++++++++++++++++Completeness in mass bins')
bin_range = [13,14.8]
nbins_x = 9
zbins = [0.2,0.5,0.8,1.0,1.15]
#zbins = [0.2,0.4,0.6,0.8,1.0,1.2]
nbins_z=len(zbins)-1
compl_m_raw = np.empty([nbins_z,nbins_x])

for i in range(0,nbins_z):
     print('-----'+str(i))
     cut1 = zbins[i]
     cut2 = zbins[i+1]
     filter1 = np.logical_and(c_merged_12.data['z_halo'] > cut1, c_merged_12.data['z_halo'] < cut2)
     c_halos_matched = c_merged_12[filter1]
     #print(c_halos_matched)
     filter2 = np.logical_and(c2.data['z'] > cut1, c2.data['z'] < cut2)
     c_halos = c2.data[filter2]
     #print(c_halos)
     h_r_halos_matched = np.histogram(c_halos_matched['log_m200c'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
     h_r_halos  = np.histogram(c_halos['log_m200c'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
     #print(h_r_halos_matched)
     #print(h_r_halos)
     compl_m_raw[i] = np.divide(h_r_halos_matched[0],h_r_halos[0],where=(h_r_halos[0]!=0))
     #print(compl_m_raw[i])

#plot
def f_completeness_param_2(log10m,log10_mc,nc):
          return np.exp(nc*np.log(10)*(log10m-log10_mc))/(1+np.exp(nc*np.log(10)*(log10m-log10_mc)))
arr_log10_mc = np.empty(nbins_z)
arr_nc = np.empty(nbins_z)
bin_x = np.empty([nbins_x])
x_bins = np.linspace(13,14.8,nbins_x+1)
labels=['0.2-0.5','0.5-0.8','0.8-1.0','1.0-1.2']
colors=['black','red','blue','purple']
#labels=['0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0','1.0-1.2']
#colors=['black','red','blue','purple','brown']
for ix in range(nbins_x):
     bin_x[ix] = 0.5 * (x_bins[ix] + x_bins[ix+1])
plt.figure()
x = np.linspace(13, 14.8, 2000)
for i in range(0,nbins_z): 
     plt.scatter(bin_x, compl_m_raw[i], label=labels[i], color=colors[i], marker= ".", s=30)
     plt.plot(bin_x, compl_m_raw[i], color=colors[i])
     popt, pcov = curve_fit(f_completeness_param_2, xdata=bin_x, ydata=compl_m_raw[i], p0=[13.5,2])
     print(popt)
     arr_log10_mc[i]=popt[0]
     arr_nc[i]=popt[1]
     f_completeness_param_2_fit = f_completeness_param_2(x, popt[0], popt[1])
     #plt.plot(x, f_completeness_param_2_fit, color=colors[i], linewidth=2.0,label="Param")
     
plt.ylim(0, 1.2)
plt.xlim(13,15)
plt.xlabel('log(m200c)')
plt.ylabel('Completeness')
plt.title(matching)
plt.legend()
plt.savefig(outpath+"completeness_mass.png", bbox_inches='tight')
plt.close()
#sys.exit()

###fit parametrizations
#define x, y and yerr
#n_z=4
#n_y=nbins_x
#n_points = n_z*n_y
#x0 = np.empty((n_points,2))
#y = np.empty(n_points)
#y_err = np.empty(n_points)
#bin_x1=np.empty([n_z])
#for ix1 in range(n_z):
#     bin_x1[ix1] = 0.5 * (zbins[ix1] + zbins[ix1+1])
#bin_x2=np.empty([n_y])
#for ix2 in range(n_y):
#     bin_x2[ix2] = 0.5 * (x_bins[ix2] + x_bins[ix1+2])
#
#for i in range(n_z):
#     for j in range(n_y):
#          print(str(i) + ' ' + str(j))
#          k = i * n_y  + j
#          print(str(k))
#          x0[k] =  np.array([bin_x1[i],bin_x2[j]])
#          y[k] = compl_m_raw[i][j]
#          y_err[k] = 0.01
#          print(str(x0[k])+' : ' +str(y[k]))
#likelihood computation
#log10_mc = 13.5
#def log_likelihood(theta, x, y, yerr):
#               c0, c1, n0, n1 = theta
#               #model = (x[:,1]/(log10_mc + c0 + c1*(1+x[:,0])))**nc/(1+(x[:,1]/(log10_mc + c0 + c1*(1+x[:,0])))**nc)
#               model = (x[:,1]/(log10_mc + c0 + c1*(1+x[:,0])))**(n0+n1*(1+x[:,0]))/(1+(x[:,1]/(log10_mc + c0 + c1*(1+x[:,0]))))**(n0+n1*(1+x[:,0]))
#               sigma2 = yerr**2
#               return -0.5*np.sum((y - model)**2/sigma2)
#
#print('Start log_likelihood minimization')
#nll = lambda *args: -log_likelihood(*args)
#initial = np.array([0,0,2.2,0.1])
##likelihood minimization minimization
#soln = minimize(nll, initial, args=(x0, y, y_err))
#c0_ml, c1_ml, n0_ml, n1_ml= soln.x
#print('likelihood results')
#print(c0_ml)
#print(c1_ml)
#print(n0_ml) 
#print(n1_ml)

#try 2-step fit
for i in range(0,nbins_z):
     print(arr_log10_mc[i])
     print(arr_nc[i])

bin_x1=np.empty([nbins_z])
for ix1 in range(nbins_z):
     bin_x1[ix1] = 0.5 * (zbins[ix1] + zbins[ix1+1]) 

def f_lin_z(z,a,b):
     return a + b*(1+z)
x = np.linspace(0.2, 1.2, 2000)

#summary plot for mc(z)
plt.figure()
plt.scatter(bin_x1, arr_log10_mc[:], label="", color= "black",marker= ".", s=30)
plt.plot(bin_x1, arr_log10_mc[:], color= "black", linewidth=1)
popt, pcov = curve_fit(f_lin_z, xdata=bin_x1, ydata=arr_log10_mc, p0=[13,0.5])
print(popt) 
a_mc=popt[0]
b_mc=popt[1]
err_a_mc=pcov[0,0]
err_b_mc=pcov[0,1]
f_lin_z_fit_mc = f_lin_z(x, a_mc, b_mc)
plt.plot(x, f_lin_z_fit_mc, color='red', linewidth=2.0,label="Fit")
plt.xlabel('z')
plt.ylabel('mc')
#plt.legend()
plt.savefig(outpath+"mc_vs_z_completeness.png", bbox_inches='tight')

#summary plot for nc(z)
plt.figure()
plt.scatter(bin_x1, arr_nc[:], label="", color= "black",marker= ".", s=30)
plt.plot(bin_x1, arr_nc[:], color= "black", linewidth=1)
popt, pcov = curve_fit(f_lin_z, xdata=bin_x1, ydata=arr_nc, p0=[13,0.5])
a_nc=popt[0]
b_nc=popt[1]
err_a_nc=pcov[0,0]
err_b_nc=pcov[0,1]
f_lin_z_fit_nc = f_lin_z(x, a_nc, b_nc)
plt.plot(x, f_lin_z_fit_nc, color='red', linewidth=2.0,label="Fit")
plt.xlabel('z')
plt.ylabel('nc')
plt.legend()
plt.savefig(outpath+"nc_vs_z_completeness.png", bbox_inches='tight')

#do completeness plots for this parametrization
plt.figure()
x = np.linspace(13, 14.8, 2000)
print(bin_x)
print(bin_x1)
for i in range(0,nbins_z):
     plt.scatter(bin_x, compl_m_raw[i], label=labels[i], color=colors[i], marker= ".", s=30)
     plt.plot(bin_x, compl_m_raw[i], color=colors[i])
     log10_mc = f_lin_z(bin_x1[i],a_mc,b_mc)
     nc = f_lin_z(bin_x1[i],a_nc,b_nc)
     f_completeness_param_2_fit = f_completeness_param_2(x, log10_mc, nc)
     plt.plot(x, f_completeness_param_2_fit, color=colors[i], linewidth=2.0,label="Param")
plt.ylim(0, 1.2)
plt.xlim(13,15)
plt.xlabel('log(m200c)')
plt.ylabel('Completeness')
plt.title(matching)
plt.legend()
plt.savefig(outpath+"completeness_mass_2_step_fit.png", bbox_inches='tight')
plt.close()
                                                   
###print results
a_nc=round(a_nc,4)
err_a_nc=round(err_a_nc,4)
b_nc=round(b_nc,4)
err_b_nc=round(err_b_nc,4)
a_mc=round(a_mc,4)
err_a_mc=round(err_a_mc,4)
b_mc=round(b_mc,4)
err_b_mc=round(err_b_mc,4)
print('Completeness parametrization: C(log10_m,z_halo) = (log10_m/log10_mc)^nc(z_halo) / ((log10_m/log10_mc)^nc(z_halo) + 1)')
print('nc(z) = a_nc + b_nc*(1+z)')
print('a_nc = ' + str(a_nc) + ' +/- ' + str(err_a_nc))
print('b_nc = ' + str(b_nc) + ' +/- ' + str(err_b_nc))
print('log10_mc(z) = a_mc + b_mc*(1+z)')
print('a_mc = ' + str(a_mc) + ' +/- ' + str(err_a_mc))
print('b_mc = ' + str(b_mc) + ' +/- ' + str(err_b_mc))

#versus richness
print('++++++++++++++++Purity in richness bins')
bin_range = [1.5,4.5]
nbins_x = 6
zbins = [0.2,0.5,0.8,1.0,1.15]
purity_m_raw = np.empty([4,nbins_x])

for i in range(0,4):
     print('-----'+str(i))
     cut1 = zbins[i]
     cut2 = zbins[i+1]
     filter1 = np.logical_and(c_merged_12.data['z_cl'] > cut1, c_merged_12.data['z_cl'] < cut2)
     c_clusters_matched = c_merged_12[filter1]
     #print(c_halos_matched)
     filter2 = np.logical_and(c1.data['z'] > cut1, c1.data['z'] < cut2)
     c_clusters = c1.data[filter2]
     #print(c_halos)
     h_r_clusters_matched = np.histogram(np.log(c_clusters_matched['richness']), bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
     h_r_clusters  = np.histogram(np.log(c_clusters['mass']), bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
     print(h_r_clusters_matched)
     print(h_r_clusters)
     purity_m_raw[i] = np.divide(h_r_clusters_matched[0],h_r_clusters[0],where=(h_r_clusters[0]!=0))
     print(purity_m_raw[i])

#plot
bin_x = np.empty([nbins_x])
x_bins = np.linspace(1.5,4.5,nbins_x+1)
arr_log10_mc = np.empty(nbins_z)
arr_nc = np.empty(nbins_z)
labels=['0.2-0.5','0.5-0.8','0.8-1.0','1.0-1.2']
colors=['black','red','blue','purple']
for ix in range(nbins_x):
     bin_x[ix] = 0.5 * (x_bins[ix] + x_bins[ix+1])
plt.figure()
for i in range(0,4): 
     plt.scatter(bin_x, purity_m_raw[i], label=labels[i], color=colors[i], marker= ".", s=30)
     plt.plot(bin_x, purity_m_raw[i], color=colors[i])
     popt, pcov = curve_fit(f_completeness_param_2, xdata=bin_x, ydata=purity_m_raw[i], p0=[1,2])
     arr_log10_mc[i]=popt[0]
     arr_nc[i]=popt[1]
     f_completeness_param_2_fit = f_completeness_param_2(x, popt[0], popt[1])
     plt.plot(x, f_completeness_param_2_fit, color=colors[i], linewidth=2.0,label="Param")
plt.ylim(0, 1.2)
plt.xlim(1.5,4.5)
plt.xlabel('ln(richness)')
plt.ylabel('Purity')
plt.title(matching)
plt.legend()
plt.savefig(outpath+"purity_richness.png", bbox_inches='tight')
plt.close()

#try 2-step fit
for i in range(0,nbins_z):
     print(arr_log10_mc[i])
     print(arr_nc[i])

bin_x1=np.empty([nbins_z])
for ix1 in range(nbins_z):
     bin_x1[ix1] = 0.5 * (zbins[ix1] + zbins[ix1+1]) 

def f_lin_z(z,a,b):
     return a + b*(1+z)
x = np.linspace(0.2, 1.2, 2000)

#summary plot for mc(z)
plt.figure()
plt.scatter(bin_x1, arr_log10_mc[:], label="", color= "black",marker= ".", s=30)
plt.plot(bin_x1, arr_log10_mc[:], color= "black", linewidth=1)
popt, pcov = curve_fit(f_lin_z, xdata=bin_x1, ydata=arr_log10_mc, p0=[13,0.5])
print(popt) 
a_mc=popt[0]
b_mc=popt[1]
err_a_mc=pcov[0,0]
err_b_mc=pcov[0,1]
f_lin_z_fit_mc = f_lin_z(x, a_mc, b_mc)
plt.plot(x, f_lin_z_fit_mc, color='red', linewidth=2.0,label="Fit")
plt.xlabel('z')
plt.ylabel('mc')
#plt.legend()
plt.savefig(outpath+"mc_vs_z_purity.png", bbox_inches='tight')

#summary plot for nc(z)
plt.figure()
plt.scatter(bin_x1, arr_nc[:], label="", color= "black",marker= ".", s=30)
plt.plot(bin_x1, arr_nc[:], color= "black", linewidth=1)
popt, pcov = curve_fit(f_lin_z, xdata=bin_x1, ydata=arr_nc, p0=[13,0.5])
a_nc=popt[0]
b_nc=popt[1]
err_a_nc=pcov[0,0]
err_b_nc=pcov[0,1]
f_lin_z_fit_nc = f_lin_z(x, a_nc, b_nc)
plt.plot(x, f_lin_z_fit_nc, color='red', linewidth=2.0,label="Fit")
plt.xlabel('z')
plt.ylabel('nc')
plt.legend()
plt.savefig(outpath+"nc_vs_z_purity.png", bbox_inches='tight')

#do purity plots for this parametrization
print(bin_x)
plt.figure()
x = np.linspace(1.7, 4, 2000)
for i in range(0,nbins_z):
     plt.scatter(bin_x, purity_m_raw[i], label=labels[i], color=colors[i], marker= ".", s=30)
     plt.plot(bin_x, purity_m_raw[i], color=colors[i])
     log10_mc = f_lin_z(bin_x1[i],a_mc,b_mc)
     nc = f_lin_z(bin_x1[i],a_nc,b_nc)
     f_completeness_param_2_fit = f_completeness_param_2(x, log10_mc, nc)
     plt.plot(x, f_completeness_param_2_fit, color=colors[i], linewidth=2.0,label="Param")
plt.ylim(0, 1.2)
plt.xlim(1.5,4.5)
plt.xlabel('ln(richness)')
plt.ylabel('Purity')
plt.title(matching)
plt.legend()
plt.savefig(outpath+"purity_mass_2_step_fit.png", bbox_inches='tight')
plt.close()
                                                   
###print results
a_nc=round(a_nc,4)
err_a_nc=round(err_a_nc,4)
b_nc=round(b_nc,4)
err_b_nc=round(err_b_nc,4)
a_mc=round(a_mc,4)
err_a_mc=round(err_a_mc,4)
b_mc=round(b_mc,4)
err_b_mc=round(err_b_mc,4)
print('Purity parametrization: P(log_r,z_cl) = (log_r/log_rc)^nc(z_cl) / ((log_r/log_rc)^nc(z_cl) + 1)')
print('nc(z) = a_nc + b_nc*(1+z)')
print('a_nc = ' + str(a_nc) + ' +/- ' + str(err_a_nc))
print('b_nc = ' + str(b_nc) + ' +/- ' + str(err_b_nc))
print('log_rc(z) = a_mc + b_rc*(1+z)')
print('a_mc = ' + str(a_mc) + ' +/- ' + str(err_a_mc))
print('b_mc = ' + str(b_mc) + ' +/- ' + str(err_b_mc))










sys.exit()
#check parametrization
#plot
bin_range = [13,14.8]
nbins_x = 9
zbins = [0.2,0.5,0.8,1.0,1.2]
bin_x = np.empty([nbins_x])
x_bins = np.linspace(13,14.8,nbins_x+1)
labels=['0.2-0.5','0.5-0.8','0.8-1.0','1.0-1.2']
colors=['black','red','blue','purple']
for ix in range(nbins_x):
          bin_x[ix] = 0.5 * (x_bins[ix] + x_bins[ix+1])
plt.figure()
for i in range(1,2):
     plt.scatter(bin_x, compl_m_raw[i], label=labels[i], color=colors[i], marker= ".", s=30)
     plt.plot(bin_x, compl_m_raw[i], color=colors[i])
plt.ylim(0, 1.2)
plt.xlim(13,15)
plt.xlabel('log(m200c)')
plt.ylabel('Completeness')
plt.title(matching)
plt.legend()
x = np.linspace(13, 14.8, 2000) 
f0 = f_completeness(x,0.65)
plt.plot(x, f0, color='red', linewidth=2.0,label="Param")
plt.savefig(outpath+"completeness_mass_check.png", bbox_inches='tight')
plt.close()

#check function alone
def f_completeness_param(log10m,log10_mc,z,c0,c1,nc):
     return np.exp(nc*np.log(10)*(log10m-(log10_mc + c0 + c1*(1+z))))/(1+np.exp(nc*np.log(10)*(log10m-(log10_mc + c0+ c1*(1+z)))))

plt.figure()
plt.ylim(0, 1.2)
plt.xlim(13,15)
plt.xlabel('log(m200c)')
plt.ylabel('Completeness')
plt.title(matching)
for i in range(0,4):
     plt.scatter(bin_x, compl_m_raw[i], label=labels[i], color=colors[i], marker= ".", s=30)
     plt.plot(bin_x, compl_m_raw[i], color=colors[i])
x = np.linspace(13, 14.8, 2000)
#f0 = f_completeness_param(x,13.4,0.65,0,0,2)
#plt.plot(x, f0, color='black', linewidth=2.0,label="Param")
#f1 = f_completeness_param(x,13.5,0.65,0,0,2)
#plt.plot(x, f1, color='blue', linewidth=2.0,label="Param")
f2 = f_completeness_param(x,13.75,1.0,0,0,2.6)
plt.plot(x, f2, color='purple', linewidth=2.0,label="Param")
plt.legend()
plt.savefig(outpath+"completeness_mass_check_functions.png", bbox_inches='tight')
plt.close()

print('DONE')
