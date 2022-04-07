#!/usr/bin/env python
# coding: utf-8

###import
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
import numpy as np
from scipy.optimize import curve_fit
#import esutil
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

###plot style
plt.rcParams['figure.figsize'] = [9.5, 6]
plt.rcParams.update({'font.size': 18})

def line(x,a,b):
             return a*x+b
def line_log(x,a,b):
             return a*np.log10(x)+b
def line_log_plot(x,a,b):
             return 10**(a*np.log10(x)+b)
def line_log_log(x,a,b):
             return 10**b*(x**a)

matching_folder = '/sps/lsst/users/tguillem/DESC/desc_may_2021/cluster_challenge/clevar_catalogs/after_matching/'

##########select case
#catalog1 = 'wazp.fits'
#catalog1 = 'redmapper.fits'
#catalog2 = 'halos.fits'
#catalog2 = 'wazp.fits'
#catalog2 = 'redmapper.fits'
catalog1 = 'c1.fits'
catalog2 = 'c2.fits'
##########

outpath = "/sps/lsst/users/tguillem/DESC/desc_may_2021/cluster_challenge/plots/"
if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)

#load c1 and c2
c1 = ClCatalog.read_full(matching_folder + catalog1)
c2 = ClCatalog.read_full(matching_folder + catalog2)
print(c1)
print(c2)

del c1['mt_cross'], c2['mt_cross']
mt1, mt2 = get_matched_pairs(c1, c2, 'cross')

#c1_cross = c1[c1['mt_cross']!=None]
#print(c1_cross)

###BUG Clevar
#del c1['mt_cross'], c2['mt_cross']
#mt = ProximityMatch()
#mt.load_matches(c1, c2, out_dir=matching_folder)
#now restrict to matched pairs
from clevar.match import get_matched_pairs
mt1, mt2 = get_matched_pairs(c1, c2, 'cross')

#count clusters
filter1 = mt1['z'] > 1.15
c_clusters_z = mt1[filter1]
print(len(mt1))
print(len(c_clusters_z))

#mass-richness plot
plt.figure()
plt.scatter(mt2['mass'], mt1['mass'], marker='.',color = 'blue', s=10, alpha=0.3, label='clusters')
#plt.scatter(np.log(mt1['mass']), np.log(mt2['mass']), marker='.',color = 'blue', s=10, alpha=0.3, label='clusters')
plt.xscale('log')
plt.yscale('log')
plt.xlim([10, 200])
plt.ylim([10, 200])
plt.xlabel('r')
#plt.ylabel('r_RM')
#plt.xlim([1, 2.5])
#plt.ylim([1, 2.5])
#plt.ylim([1.0e13, 2.0e15])
#plt.ylim([10**13, 300])
#plt.ylim([10, 300])
plt.xlabel('r_RM')
plt.ylabel('r_WaZP')
#plt.xlabel('r_WaZP')
#plt.ylabel('mass')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
#a, b = np.polyfit(np.log10(mt1['mass']), np.log10(mt1['mass']), 1)
#x = np.linspace(1, 300, 1000)
#line1 = line(x, a, b)
#line1 = line_log_log(x, a, b)
#plt.plot(x, line1, color='red', linewidth=0.5)
#x = np.logspace(10, 300, 1000)
#norm = m*100+b/(10**13)
#print(a)
#print(b)
#fit 2 with log
popt, pcov = curve_fit(line, xdata=np.log10(mt1['mass']), ydata=np.log10(mt2['mass']), p0=[0.9, 0.01])
print(popt)
x = np.linspace(3, 300, 1000)
line2 = line_log_plot(x, popt[0], popt[1])
plt.plot(x, line2, color='red', linewidth=0.5)
plt.savefig(outpath+'mass_richness.png', bbox_inches='tight')
plt.close()
#sys.exit()

###test log-log plot
plt.figure()
plt.scatter(np.log10(mt1['mass']), np.log10(mt2['mass']), marker='.',color = 'blue', s=10, alpha=0.3, label='clusters')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0.5, 2.5])
plt.ylim([0.5, 2.5])
plt.xlabel('log(r_WaZP)')
plt.ylabel('log(r_RM)')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
#fit without log
#a, b = np.polyfit(mt1['mass'], mt2['mass'], 1)
#x = np.linspace(10, 300, 1000)
x = np.linspace(0.5, 2.5, 1000)
line1 = line(x, popt[0], popt[1])
#line1 = line_log_log(x, a, b)
plt.plot(x, line1, color='red', linewidth=0.5)
#x = np.logspace(10, 300, 1000)
#norm = m*100+b/(10**13)
#print(a)
#print(b)
plt.savefig(outpath+'mass_richness_log_log.png', bbox_inches='tight')
plt.close()
#sys.exit()

###clevar test
plt.figure()
info = scaling.mass_density_metrics(
        c2[c2['mass']>20], c1, 'cross', ax_rotation=45,
        add_fit=False, fit_bins1=8)
plt.savefig(outpath+'mass_richness_clevar.png', bbox_inches='tight')
plt.close()

#create unmatched catalogs
wazp_unmatched = c1[c1['mt_cross']==None]
redmapper_unmatched = c2[c2['mt_cross']==None]
print(wazp_unmatched)
print(redmapper_unmatched)
print(redmapper_unmatched['id'])
redmapper_unmatched.write(outpath + 'redmapper_unmatched.fits', overwrite=True)

#redshift
plt.figure()
plt.scatter(mt2['z'], mt1['z'], marker='.',color = 'blue', s=10, alpha=0.3, label='clusters')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0,1.5])
plt.ylim([0,1.5])
plt.xlabel('z_RM')
plt.ylabel('z_WaZP')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'redshift_cluster_halo.png', bbox_inches='tight')
plt.close()
#sys.exit()

#check redshift
plt.figure()
bin_range = [0,1.5]
nbins = 15
#plt.hist(galaxy_data_2['redshift'], range=bin_range, bins=nbins, label='m_z*', histtype='step', color = 'blue')
plt.hist(c1['z'], range=bin_range, bins=nbins, label='all W', histtype='step', color = 'blue')
plt.hist(c2['z'], range=bin_range, bins=nbins, label='all R', histtype='step', color = 'red') 
plt.xlabel("z");
plt.ylabel("clusters")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
#plt.title('')
plt.legend(title = '', loc='upper right')
plt.savefig(outpath+'clusters_vs_z.png')
plt.close() 

#check richness
plt.figure()
bin_range = [0,100]
nbins = 20
#plt.hist(galaxy_data_2['redshift'], range=bin_range, bins=nbins, label='m_z*', histtype='step', color = 'blue')
plt.hist(c1['mass'], range=bin_range, bins=nbins, label='all W', histtype='step', color = 'blue')
plt.hist(c2['mass'], range=bin_range, bins=nbins, label='all R', histtype='step', color = 'red') 
plt.xlabel("z");
plt.ylabel("clusters")
#plt.yscale('log')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
#plt.title('')
plt.legend(title = '', loc='upper right')
plt.savefig(outpath+'clusters_vs_richness.png')
plt.close()

#check redshift
plt.figure()
bin_range = [0,1.5]
nbins = 15
#plt.hist(galaxy_data_2['redshift'], range=bin_range, bins=nbins, label='m_z*', histtype='step', color = 'blue')
plt.hist(mt1['z'], range=bin_range, bins=nbins, label='cross W', histtype='step', color = 'blue')
plt.hist(mt2['z'], range=bin_range, bins=nbins, label='cross R', histtype='step', color = 'red') 
plt.xlabel("z");
plt.ylabel("clusters")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
#plt.title('')
plt.legend(title = '', loc='upper right')
plt.savefig(outpath+'clusters_vs_z_cross.png')
plt.close() 

#check richness
plt.figure()
bin_range = [0,100]
nbins = 20
#plt.hist(galaxy_data_2['redshift'], range=bin_range, bins=nbins, label='m_z*', histtype='step', color = 'blue')
plt.hist(mt1['mass'], range=bin_range, bins=nbins, label='cross W', histtype='step', color = 'blue')
plt.hist(mt2['mass'], range=bin_range, bins=nbins, label='cross R', histtype='step', color = 'red') 
plt.xlabel("r");
plt.ylabel("clusters")
#plt.yscale('log')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
#plt.title('')
plt.legend(title = '', loc='upper right')
plt.savefig(outpath+'clusters_vs_richness_cross.png')
plt.close()

#check redshift
plt.figure()
bin_range = [0,1.5]
nbins = 15
#plt.hist(galaxy_data_2['redshift'], range=bin_range, bins=nbins, label='m_z*', histtype='step', color = 'blue')
plt.hist(wazp_unmatched['z'], range=bin_range, bins=nbins, label='excl W', histtype='step', color = 'blue')
plt.hist(redmapper_unmatched['z'], range=bin_range, bins=nbins, label='excl R', histtype='step', color = 'red') 
plt.xlabel("z");
plt.ylabel("clusters")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
#plt.title('')
plt.legend(title = '', loc='upper right')
plt.savefig(outpath+'clusters_vs_z_excl.png')
plt.close() 

#check richness
plt.figure()
bin_range = [0,100]
nbins = 20
#plt.hist(galaxy_data_2['redshift'], range=bin_range, bins=nbins, label='m_z*', histtype='step', color = 'blue')
plt.hist(wazp_unmatched['mass'], range=bin_range, bins=nbins, label='excl W', histtype='step', color = 'blue')
plt.hist(redmapper_unmatched['mass'], range=bin_range, bins=nbins, label='excl R', histtype='step', color = 'red') 
plt.xlabel("r");
plt.ylabel("clusters")
#plt.yscale('log')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
#plt.title('')
plt.legend(title = '', loc='upper right')
plt.savefig(outpath+'clusters_vs_richness_excl.png')

#sys.exit()

##########################create a merged catalog
from clevar.match import output_matched_catalog

#output_matched_catalog(matching_folder+'wazp.fits', matching_folder+'redmapper.fits',matching_folder+'output_catalog.fits', c1, c2, matching_type='cross', overwrite=True)
#output_matched_catalog(matching_folder+'wazp.fits', matching_folder+'halos.fits',matching_folder+'output_catalog.fits', c1, c2, matching_type='cross', overwrite=True)
output_matched_catalog(matching_folder+catalog1, matching_folder+catalog2,matching_folder+'output_catalog.fits', c1, c2, matching_type='cross', overwrite=True)
#test multi
#output_matched_catalog(matching_folder+catalog1, matching_folder+catalog2,matching_folder+'output_catalog.fits', c1, c2, matching_type='multi_other', overwrite=True)
#cat_test = Table.read(matching_folder+'output_catalog.fits') 
#print(cat_test)
c_merged = ClCatalog.read(matching_folder+'output_catalog.fits', 'merged',  cat1_z='cat1_z', cat1_mass = 'cat1_mass', cat2_z='cat2_z', cat2_mass = 'cat2_mass')
#c_merged = ClCatalog.read(matching_folder+'output_catalog.fits', 'merged', cat1_ra='cat1_SkyCoord[0]', cat1_dec='cat1_dec', cat1_z='cat1_z', cat1_mass = 'cat1_mass', cat2_ra='cat2_ra', cat2_dec='cat2_dec', cat2_z='cat2_z', cat2_mass = 'cat2_mass')

zbins = [0,0.5,0.75,1.0,1.25,1.5]
ybins = 10**np.linspace(13, 15, 20)
xbins = [25,30,35,40,50,60,100]

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
print('mass 2D')

#redshift from merged catalog
plt.figure()
plt.scatter(c_merged.data['cat1_z'], c_merged.data['cat2_z'], marker='.',color = 'blue', s=10, alpha=0.3, label='clusters')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0,1.5])
plt.ylim([0,1.5])
plt.xlabel('z_cat1')
plt.ylabel('z_cat2')
x = np.linspace(0, 1.5, 1000)
line1 = line(x, 1, 0)
plt.plot(x, line1, color='red', linewidth=0.5)
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'redshift_cluster_halo_fix.png', bbox_inches='tight')
plt.close()
#debug delta_z
#merged_delta_z = c_merged.data[(c_merged.data['cat1_z']-c_merged.data['cat2_z'])>0.02*(1+c_merged.data['cat2_z'])]
#print(merged_delta_z['cat1_z','cat2_z'])
#print(merged_delta_z['cat1_z','cat2_z',math.fabs(c_merged.data['cat1_z']-c_merged.data['cat2_z']),0.05*(1+c_merged.data['cat2_z'])])

###study the unmatched clusters of wazp
#wazp_unmatched = Table(c1['mass'])#,c1['mt_cross'])#,filters=[c1.data['mt_cross']==None])
#print(wazp_unmatched)
wazp_unmatched = c1[c1.data['mt_cross']==None]
#print(wazp_unmatched)
redmapper_unmatched = c2[c2.data['mt_cross']==None]
print(redmapper_unmatched['id'])

#for sky area normalization
entries = len(c2.data['mass'])
print(entries)
dm = np.empty([entries])
for ix in range(entries):
     #dm[ix] = 0.13
     dm[ix] = 0.4 #54/135
plt.figure()
bin_range = [0,1.5]
nbins = 15
#bin_range = [1.100,1.150]
#nbins = 50
plt.hist(c1.data['z'], range=bin_range, bins=nbins, label='all', histtype='step', color = 'red')
plt.hist(mt1['z'], range=bin_range, bins=nbins, label='matched', histtype='step', color = 'black')
plt.hist(c2.data['z'], range=bin_range, bins=nbins, label='RM', histtype='step', color = 'blue', weights=dm)# density=True) 
plt.xlabel("z");
plt.ylabel("clusters")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.title('wazp cosmoDC2: matching to wazp')
plt.legend(title = '', loc='upper left')
plt.savefig(outpath+'clusters_vs_z.png')
plt.close()

#plot richness
plt.figure() 
bin_range = [20,100]
nbins = 16
plt.hist(c1.data['mass'], range=bin_range, bins=nbins, label='all', histtype='step', color = 'red')
plt.hist(mt1['mass'], range=bin_range, bins=nbins, label='matched', histtype='step', color = 'black')
plt.hist(c2.data['mass'], range=bin_range, bins=nbins, label='RM', histtype='step', color = 'blue', weights=dm)#, density=True)
plt.xlabel("NGALS");
plt.ylabel("clusters")
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.title('wazp cosmoDC2: matching to redmapper')
plt.legend(title = '', loc='upper right')
plt.savefig(outpath+'clusters_vs_mass.png')
plt.close()

###completeness
print('++++++++++++++++Completeness')
#print(c2)
#c1.data, mt1, cluster_unmatched
#c2.data, mt2, halos_unmatched
cluster_unmatched = c1[c1.data['mt_cross']==None]
#print(wazp_unmatched)
halos_unmatched = c2[c2.data['mt_cross']==None]
#halos_matched_multi = c2[c2.data['cat2_mt_multi_self']!=None] 
#versus NGALS
#h_NGALS_wazp_matched = np.histogram(mt1['mass'], bins=nbins, range=bin_range, normed=None, weights=None, density=None)
#h_halos = np.histogram(c2.data['mass'], bins=nbins, range=bin_range, normed=None, weights=None, density=None)
#completeness_mass = np.empty([entries])
#for ix in range(nbins):
#     print(h_NGALS_wazp_matched[0][ix])
#     print(h_halos[0][ix])

#versus z
bin_range = [0.2,1.1]
nbins_x = 9 
#nbins_y = 10
x_bins = np.linspace(0.2,1.1,nbins_x+1)
mbins = [10**13,10**13.5,10**14,10**14.5,10**15]
compl_z_raw = np.empty([4,nbins_x])
compl_z = np.empty([4,nbins_x])

for i in range(0,4):
    #print(i)
    cut1 = mbins[i]
    cut2 = mbins[i+1]
    filter1 = np.logical_and(mt2['mass'] > cut1, mt2['mass'] < cut2)
    c_halos_matched = mt2[filter1]
    filter2 = np.logical_and(c2.data['mass'] > cut1, c2.data['mass'] < cut2)
    c_halos = c2.data[filter2]
    #filter3 = np.logical_and(halos_matched_multi['mass'] > cut1, halos_matched_multi['mass'] < cut2)
    #c_halos_matched_multi = halos_matched_multi[filter3]
    h_z_halos_matched = np.histogram(c_halos_matched['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    #multi
    #h_z_halos_matched = np.histogram(c_halos_matched_multi['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None) 
    h_z_halos = np.histogram(c_halos['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    n1_halo_min = 10
    number1_of_match = h_z_halos_matched[0]
    number1_of_halo = h_z_halos[0]
    #print(number1_of_match)
    #print(number1_of_halo)
    compl_z_raw[i] = (number1_of_match/number1_of_halo)
    #compl_z[i] = np.ma.masked_where(number1_of_halo<n1_halo_min, compl_z_raw)
    print(compl_z_raw[i])
    #print(compl_z[i])
    
print('loop ok')
#h_z_halos_matched = np.histogram(mt2['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
#h_z_halos = np.histogram(c2.data['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
#completeness_z = np.empty([nbins])
#for ix in range(nbins_x):
#     print(h_z_cluster_matched[0][ix])
#     print(h_z_halos[0][ix])
#     completeness_z[ix] = h_z_cluster_matched[0][ix] / h_z_halos[0][ix]
#     print("completeness = " + str(completeness_z[ix]))

#n1_halo_min = 10
#number1_of_match = h_z_halos_matched[0]
#number1_of_halo = h_z_halos[0]
#compl_z_raw = (number1_of_match/number1_of_halo)
#compl_z = np.ma.masked_where(number1_of_halo<n1_halo_min, compl_z_raw)
#print(compl_z)
plt.figure()
bin_x = np.empty([nbins_x])
print(nbins_x)
for ix in range(nbins_x):
    bin_x[ix] = 0.5 * (x_bins[ix] + x_bins[ix+1])
plt.scatter(bin_x, compl_z_raw[0], label="13-13.5", color= "black",marker= ".", s=30)
plt.plot(bin_x, compl_z_raw[0], color= "black")
plt.scatter(bin_x, compl_z_raw[1], label="13.5-14", color= "red",marker= ".", s=30)
plt.plot(bin_x, compl_z_raw[1], color= "red")
plt.scatter(bin_x, compl_z_raw[2], label="14-14.5", color= "blue",marker= ".", s=30)
plt.plot(bin_x, compl_z_raw[2], color= "blue")
plt.scatter(bin_x, compl_z_raw[3], label="14.5-15", color= "purple",marker= ".", s=30)
plt.plot(bin_x, compl_z_raw[3], color= "purple")
plt.ylim(0,1.2)
plt.xlim(0.2,1.6)
plt.xlabel('z')
plt.ylabel('Completeness')
plt.title(algo)
plt.legend()
#plt.legend(loc = (1.05,0.7), prop = {'size': 15})
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+"completeness_z.png", bbox_inches='tight')
plt.close()
#sys.exit()

#versus richness
bin_range = [0,100]
nbins_x = 10
mbins = [10**13,10**13.5,10**14,10**14.5,10**15]
compl_r_raw = np.empty([4,10])
compl_r = np.empty([4,10])

for i in range(0,4):
    print('------ ' + str(i))
    cut1 = mbins[i]
    cut2 = mbins[i+1]
    filter1 = np.logical_and(c_merged.data['cat2_mass'] > cut1, c_merged.data['cat2_mass'] < cut2)
    c_halos_matched = c_merged[filter1]
    filter2 = np.logical_and(c2.data['mass'] > cut1, c2.data['mass'] < cut2)
    c_halos = c2.data[filter2]
    #filter3 = np.logical_and(halos_matched_multi['mass'] > cut1, halos_matched_multi['mass'] < cut2)
    #c_halos_matched_multi = halos_matched_multi[filter3]
    h_r_halos_matched = np.histogram(c_halos_matched['cat1_mass'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    #multi
    #h_z_halos_matched = np.histogram(c_halos_matched_multi['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None) 
    #h_r_halos = np.histogram(c_halos['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    n1_halo_min = 10
    number1_of_match = h_r_halos_matched[0]
    number1_of_halo = len(c_halos)
    print(number1_of_match)
    print(number1_of_halo)
    compl_r_raw[i] = (number1_of_match/number1_of_halo)
    #compl_z[i] = np.ma.masked_where(number1_of_halo<n1_halo_min, compl_z_raw)
    print(compl_r_raw[i])
    #print(compl_z[i])
    
print('loop richnes ok')
#h_z_halos_matched = np.histogram(mt2['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
#h_z_halos = np.histogram(c2.data['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
#completeness_z = np.empty([nbins])
#for ix in range(nbins_x):
#     print(h_z_cluster_matched[0][ix])
#     print(h_z_halos[0][ix])
#     completeness_z[ix] = h_z_cluster_matched[0][ix] / h_z_halos[0][ix]
#     print("completeness = " + str(completeness_z[ix]))

#n1_halo_min = 10
#number1_of_match = h_z_halos_matched[0]
#number1_of_halo = h_z_halos[0]
#compl_z_raw = (number1_of_match/number1_of_halo)
#compl_z = np.ma.masked_where(number1_of_halo<n1_halo_min, compl_z_raw)
#print(compl_z)
#compute compl_r_raw_incl
plt.figure()
bin_x = np.empty([nbins_x])
x_bins = np.linspace(0,100,10)
for ix in range(nbins_x-1):
    bin_x[ix] = 0.5 * (x_bins[ix] + x_bins[ix+1])
plt.scatter(bin_x, compl_r_raw[0], label="13-13.5", color= "black",marker= ".", s=30)
plt.plot(bin_x, compl_r_raw[0], color= "black")
plt.scatter(bin_x, compl_r_raw[1], label="13.5-14", color= "red",marker= ".", s=30)
plt.plot(bin_x, compl_r_raw[1], color= "red")
plt.scatter(bin_x, compl_r_raw[2], label="14-14.5", color= "blue",marker= ".", s=30)
plt.plot(bin_x, compl_r_raw[2], color= "blue")
plt.scatter(bin_x, compl_r_raw[3], label="14.5-15", color= "purple",marker= ".", s=30)
plt.plot(bin_x, compl_r_raw[3], color= "purple")
plt.ylim(0, 1.2)
plt.xlim(0,100)
plt.xlabel('r')
plt.ylabel('Completeness')
plt.title(algo)
plt.legend()
plt.savefig(outpath+"completeness_r.png", bbox_inches='tight')
plt.close()
#sys.exit()

#2D : z-M plane
#bin_range = [0.2,1.1]
nbins_x = 9 
nbins_y = 13
x_bins = np.linspace(0.2,1.1,10)
#nbins_x = 16
##nbins_y = 13
##x_bins = np.linspace(0,1.5,16)
#y_bins = np.logspace(10**13,10**15,10)
y_bins = np.linspace(13,15.4,13)
#print(x_bins)
#print(y_bins)
h2_z_halos_matched = np.histogram2d(mt2['z'], np.log10(mt2['mass']), bins=[nbins_x,nbins_y], range=[[x_bins[0],x_bins[nbins_x]],[13,15]], normed=None, weights=None, density=None)
h2_z_halos = np.histogram2d(c2.data['z'], np.log10(c2.data['mass']), bins=[nbins_x,nbins_y], range=[[x_bins[0],x_bins[nbins_x]],[13,15]], normed=None, weights=None, density=None)

n_halo_min = 10
number_of_match = h2_z_halos_matched[0]
number_of_halo = h2_z_halos[0]
compl_2d_raw = (number_of_match/number_of_halo)
compl_2d = np.ma.masked_where(number_of_halo<n_halo_min, compl_2d_raw)
#print(compl_2d)
#print("2D")

#NOT USED FOR NOW: bin-by-bin computation of the completeness
completeness_z_halo_mass = [[0]*(nbins_y)]*(nbins_x)
for ix in range(nbins_x):
    #print(ix)
    for iy in range(nbins_y):
        #print("loop iy")
        #print(iy)
        #print((h2_z_halos_matched[0][ix])[iy])
        #print((h2_z_halos[0][ix])[iy])
        inclusive = (h2_z_halos[0][ix])[iy]
        matched = (h2_z_halos_matched[0][ix])[iy]
        completeness_z_halo_mass[ix][iy]=0
        if(inclusive!=0):
            completeness_z_halo_mass[ix][iy] = matched/inclusive
            #print("completeness = " + str(completeness_z_halo_mass[ix][iy]))
#NOT USED FOR NOW

#2D plot
fig, ax = plt.subplots()
x, y = np.meshgrid(x_bins, y_bins)
#print(x)
#print(y)
c = ax.pcolormesh(x, y, compl_2d.T, cmap='jet', vmin=0, vmax=1)
ax.set_title('Completeness 2D: ' + algo)
ax.set_xlim(0.2,1.1)
ax.set_ylim(13.0,15.4)
ax.set_xlabel('z')
ax.set_ylabel('log(halo_mass)') 
fig.colorbar(c, ax=ax)
fig.savefig(outpath+'completeness_2D.png')
#sys.exit()

#density
#plt.figure()
bin_range = [13,15]
nbins = 4
plt.hist(np.log10(c_merged.data['cat2_mass']), range=bin_range, bins=nbins, label='DC2', histtype='step', color = 'black')
plt.xlabel("halo_mass");
plt.ylabel("halos")
plt.yscale('log')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.title('cosmoDC2')
plt.legend(title = '', loc='upper right')
plt.savefig(outpath+'halos_vs_mass.png')
plt.close() 
#sys.exit()

print('++++++++++++++++Purity')     
#versus z
bin_range = [0,1.5]
nbins_x = 15 
nbins_y = 10
x_bins = np.linspace(0,1.5,16)
#mbins = [10**13,10**13.5,10**14,10**14.5,10**15]
#wazp
mbins = [25,40,70,200]
#redmapper
#mbins = [25,30,60,200]
purity_z_raw = np.empty([3,15])
purity_z = np.empty([3,15])

scaling_for_RM = 2.37
for i in range(0,3):
    print(i)
    cut1 = mbins[i]
    cut2 = mbins[i+1]
    filter1 = np.logical_and(mt1['mass'] > cut1, mt1['mass'] < cut2)
    c_clusters_matched = mt1[filter1]
    filter2 = np.logical_and(c1.data['mass'] > cut1, c1.data['mass'] < cut2)
    c_clusters = c1.data[filter2]
    h_z_clusters_matched = np.histogram(c_clusters_matched['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    h_z_clusters = np.histogram(c_clusters['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
    n1_halo_min = 10
    number1_of_match = h_z_clusters_matched[0]
    number1_of_halo = h_z_clusters[0]#/1.9
    print(number1_of_match)
    print(number1_of_halo)
    purity_z_raw[i] = (number1_of_match/number1_of_halo)
    #purity_z[i] = np.ma.masked_where(number1_of_halo<n1_halo_min, purity_z_raw)
    print(purity_z_raw[i])
    #print(purity_z[i])
    
print('loop ok')
#h_z_halos_matched = np.histogram(mt2['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
#h_z_halos = np.histogram(c2.data['z'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
#completeness_z = np.empty([nbins])
#for ix in range(nbins_x):
#     print(h_z_cluster_matched[0][ix])
#     print(h_z_halos[0][ix])
#     completeness_z[ix] = h_z_cluster_matched[0][ix] / h_z_halos[0][ix]
#     print("completeness = " + str(completeness_z[ix]))

#n1_halo_min = 10
#number1_of_match = h_z_halos_matched[0]
#number1_of_halo = h_z_halos[0]
#purity_z_raw = (number1_of_match/number1_of_halo)
#purity_z = np.ma.masked_where(number1_of_halo<n1_halo_min, purity_z_raw)
#print(purity_z)
plt.figure()
bin_x = np.empty([nbins_x])
for ix in range(nbins_x):
    bin_x[ix] = 0.5 * (x_bins[ix] + x_bins[ix+1])
plt.scatter(bin_x, purity_z_raw[0], label=str(mbins[0])+'-'+str(mbins[1]), color= "black",marker= ".", s=30)
plt.plot(bin_x, purity_z_raw[0], color= "black")
plt.scatter(bin_x, purity_z_raw[1], label=str(mbins[1])+'-'+str(mbins[2]), color= "red",marker= ".", s=30)
plt.plot(bin_x, purity_z_raw[1], color= "red")
plt.scatter(bin_x, purity_z_raw[2], label=str(mbins[2])+'-'+str(mbins[3]), color= "blue",marker= ".", s=30)
plt.plot(bin_x, purity_z_raw[2], color= "blue")
plt.ylim(0, 1.2)
plt.xlim(0.2,1.2)
plt.xlabel('z')
plt.ylabel('purity')
plt.title(algo)
plt.legend()
plt.savefig(outpath+"purity_z.png", bbox_inches='tight')
plt.close()



sys.exit()
