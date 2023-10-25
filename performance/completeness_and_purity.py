    
import numpy as np
import sys
import os
from astropy.table import Table
from astropy.io import fits
from astropy.io import ascii
import GCRCatalogs

from clevar.catalog import ClCatalog
from clevar.match_metrics.recovery import ClCatalogFuncs as r_cf
from clevar.match_metrics import recovery
from clevar.match import output_matched_catalog

import matplotlib.pyplot as plt
import matplotlib as mpl

from tqdm import tqdm



cats = [sys.argv[1].split('_'), sys.argv[2].split('_')]
matching_method = sys.argv[3]

(first, last) = np.argsort([c[0] for c in cats])

    
area = [440 if 'cosmoDC2' in c else 300 if 'DC2' in c else None for c in (sys.argv[1], sys.argv[2])]


inpath = '/sps/lsst/groups/clusters/cluster_comparison_project/after_matching/'
inpath += f"{'.'.join(cats[first][:-1])}_{'.'.join(cats[last][:-1])}/"		## EXAMPLE INPUT: cosmoDC2_v1 wazp_cosmoDC2.fzb_v0 proximity
inpath += f"{cats[first][-1]}_{cats[last][-1]}/"				## RESULTING INPATH: /cosmoDC2_wazp.cosmoDC2.fzb/v1_v0/
if matching_method == 'proximity' :
	inpath += 'proximity_matching_deltaz_0.05_matchradius_1mpc/'
else :
	inpath += 'member_matching_fshare_0.0/'

outpath = '/sps/lsst/users/namourou/web/desc/cluster_comparison_project/'
outpath += f"{'.'.join(cats[first][:-1])}_{'.'.join(cats[last][:-1])}/"
outpath += f"{cats[first][-1]}_{cats[last][-1]}/"
if matching_method == 'proximity' :
	outpath += 'proximity_matching_deltaz_0.05_matchradius_1mpc/'
else :
	outpath += 'member_matching_fshare_0.0/'

## READ IN MATCHED CATALOGS
cl = []
for c in cats :
	cl.append(ClCatalog.read_full(inpath + '.'.join(c) + '.fits'))


## IF CATALOGS EXIST, MAKE CORRESPONDING PLOT DIRECTORIES
outpaths = [outpath+'png_plots/', outpath+'pdf_plots/']
for path in outpaths :
	if not os.path.exists(path) :
		os.makedirs(path)
outpath+='png_plots/'
#cl[0].sort('mt_cross')
#cl[1].sort('id_cl')

#matches = [cl[0]['mt_cross'] != '', cl[1]['mt_cross'] != '']

linestyle = ('-','--') if cats[0][0]==cats[1][0] else ('-','-')

figx=10
figy=7

## CREATE A MERGED CATALOG FOR THE CROSS-MATCHED PAIRS
output_matched_catalog(inpath + '.'.join(cats[0]) + '.fits', inpath + '.'.join(cats[1]) + '.fits', inpath + 'output_catalog_12.fits', cl[0], cl[1], matching_type='cross', overwrite=True)
c_merged_12 = ClCatalog.read(inpath+'output_catalog_12.fits', name = 'merged',  full = True)
#c_merged_12 = c_merged_12[c_merged_12.data['z_cl']<1.15]
 

## COMPLETENESS VERSUS Z

if True:
    print('--- Completeness VS redshift ---')
    nbins_x = 10
    bin1 = np.linspace(0.0, 1.8, nbins_x) #For AMICO
    bin2 = [10**13,10**13.5,10**14,10**14.3,10**15]
    labels=['$10^{13}-10^{13.5}$','$10^{13.5}-10^{14}$', '$10^{14}-10^{14.3}$','$10^{14.3}-10^{15}$']
    plt.xlim(0,2.0) #For AMICO
    nbins_x -= 1
    bin_range = [min(bin1), max(bin1)]
    compl = np.empty([len(bin2),nbins_x])
    bin_x = np.empty([nbins_x])
    for ix in range(nbins_x):
        bin_x[ix] = 0.5 * (bin1[ix] + bin1[ix+1])
    for i in range(0,len(bin2)-1):
        cut1 = bin2[i]
        cut2 = bin2[i+1]
        filter1 = np.logical_and(c_merged_12.data["cat1_mass"] > cut1, c_merged_12.data["cat1_mass"] < cut2)
        filter2 = np.logical_and(cl[0].data["mass"] > cut1, cl[0].data["mass"] < cut2)
        c_halos_matched = c_merged_12[filter1]
        c_halos = cl[0][filter2]
        h_r_halos_matched = np.histogram(c_halos_matched["cat1_z_cl"], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
        h_r_halos  = np.histogram(c_halos["z_cl"], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
        compl[i] = np.divide(h_r_halos_matched[0], h_r_halos[0], where=(h_r_halos[0]!=0))
        for j in range(len(compl[i])):
            if h_r_halos_matched[0][j]<10 or h_r_halos[0][j]<10:
                compl[i][j] = np.nan
        plt.ylim(0, 1.2)
        plt.xlabel('$z$', fontsize = 13)
        plt.ylabel('Completeness', fontsize = 13)
        plt.plot(bin_x, compl[i], marker = '+', label = labels[i])

    plt.title('AMICO-cosmoDC2')
    plt.legend()
    plt.savefig(outpath + 'completeness_vs_z_mass_binned' + '.png', format='png')
    plt.close()


## COMPLETENESS VERSUS MASS

if True : 
    print('--- Completeness VS mass ---')
    nbins_x = 10
    bin1 = np.linspace(13, 15, nbins_x)
    bin2 = [.2,.5,.8,1,1.2,1.5,1.8] #For AMICO
    labels=['0.2-0.5','0.5-0.8','0.8-1.0','1.0-1.2','1.2-1.5','1.5-1.8']
    plt.xlim(13,15)

    nbins_x -= 1
    bin_range = [min(bin1), max(bin1)]
    compl = np.empty([len(bin2),nbins_x])
    bin_x = np.empty([nbins_x])
    for ix in range(nbins_x):
        bin_x[ix] = 0.5 * (bin1[ix] + bin1[ix+1])
    for i in range(0,len(bin2)-1):
        cut1 = bin2[i]
        cut2 = bin2[i+1]
        filter1 = np.logical_and(c_merged_12.data["cat1_z_cl"] > cut1, c_merged_12.data["cat1_z_cl"] < cut2)
        filter2 = np.logical_and(cl[0]["z_cl"] > cut1, cl[0]["z_cl"] < cut2)
        c_halos_matched = c_merged_12[filter1]
        c_halos = cl[0][filter2]
        h_r_halos_matched = np.histogram(c_halos_matched["cat1_log10_mass"], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
        h_r_halos  = np.histogram(c_halos["log10_mass"], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
        compl[i] = np.divide(h_r_halos_matched[0], h_r_halos[0], where=(h_r_halos[0]!=0))
        for j in range(len(compl[i])):
            if h_r_halos_matched[0][j]<10 or h_r_halos[0][j]<10:
                compl[i][j] = np.nan
        plt.ylim(0, 1.2)
        plt.xlabel('mass', fontsize = 13)
        plt.ylabel('Completeness', fontsize = 13)
        plt.plot(bin_x, compl[i], marker = '+', label = labels[i])

    plt.title('AMICO-cosmoDC2', fontsize = 13)
    plt.legend()
    plt.savefig(outpath + 'completeness_vs_mass_z_binned' + '.png')
    plt.close()
                
        
## PURITY VERSUS RICHNESS

if True : 
    print('--- Purity VS richness ---')
    bin_range = [0,100]
    nbins_x = 10
    zbins = [0,0.5,0.8,1.0,1.2,1.8] #Amico
    purity_z_raw = np.empty([len(zbins),nbins_x])

    bin_x = np.empty([nbins_x])
    x_bins = np.linspace(0,100,nbins_x+1)
    labels=['0-0.5','0.5-0.8','0.8-1.0','1.0-1.2', '1.2-1.8']

    for ix in range(nbins_x):
         bin_x[ix] = 0.5 * (x_bins[ix] + x_bins[ix+1])

    for i in range(0,len(zbins)-1):
        cut1 = zbins[i]
        cut2 = zbins[i+1]
        filter1 = np.logical_and(c_merged_12.data['cat2_z_cl'] > cut1, c_merged_12.data['cat2_z_cl'] < cut2)
        c_clusters_matched = c_merged_12[filter1]
        #print(c_clusters_matched)
        filter2 = np.logical_and(cl[1]['z_cl'] > cut1, cl[1]['z_cl'] < cut2)
        c_clusters = cl[1][filter2]
        #print(c_clusters)
        h_r_clusters_matched = np.histogram(c_clusters_matched['cat2_mass'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
        h_r_clusters  = np.histogram(c_clusters['mass'], bins=nbins_x, range=bin_range, normed=None, weights=None, density=None)
        #print(h_r_clusters_matched)
        #print(h_r_clusters)
        purity_z_raw[i] = np.divide(h_r_clusters_matched[0],h_r_clusters[0],where=(h_r_clusters[0]!=0))
        for j in range(len(purity_z_raw[i])):
            if h_r_clusters_matched[0][j]<10 or h_r_clusters[0][j]<10:
                purity_z_raw[i][j] = np.nan
        plt.scatter(bin_x, purity_z_raw[i], label=labels[i], marker= ".", s=30)
        plt.plot(bin_x, purity_z_raw[i])
    plt.xlabel('$\lambda^*$', fontsize = 13)
    plt.ylabel('Purity', fontsize = 13)
    plt.legend()
    plt.ylim(0,1.2)
    plt.xlim(0,110) 
    plt.title('AMICO-cosmoDC2', fontsize = 13)
    plt.savefig(outpath + 'purity_vs_richness_z_binned' + '.png', bbox_inches='tight', format='png')
    plt.close()

## 2D COMPLETENESS

if True :
    print('--- 2D Completeness ---')
    nbins_x = 18
    nbins_y = 21
    x_bins = np.linspace(0,1.8,nbins_x)
    y_bins = np.logspace(13,15,nbins_y)
    xbin_range = [min(x_bins), max(x_bins)]
    ybin_range = [min(y_bins), max(y_bins)]
    h2_z_halos_matched = np.histogram2d(c_merged_12['cat1_z_cl'], c_merged_12['cat1_mass'], bins=(x_bins, y_bins), 
                                        range=[xbin_range, ybin_range], normed=None, weights=None, density=None)
    h2_z_halos = np.histogram2d(cl[0]['z_cl'], cl[0]['mass'], bins=(x_bins, y_bins),
                                range=[xbin_range, ybin_range], normed=None, weights=None, density=None)
    number_of_match = h2_z_halos_matched[0]
    number_of_halo = h2_z_halos[0]
    compl_2d = np.divide(number_of_match, number_of_halo, where=(number_of_halo!=0))
    fig, ax = plt.subplots(figsize =(12,6))
    x, y = np.meshgrid(x_bins, y_bins)
    #to avoid seeing 0 values
    compl_2d[compl_2d==0] = np.nan
    c = ax.pcolormesh(x, y, compl_2d.T, cmap='jet', vmin=0, vmax=1)
    ax.set_xlim(0,2.0)
    ax.set_ylim(10**13,10**15)
    ax.set_xlabel('z', fontsize = 13)
    ax.set_ylabel('halo_mass', fontsize = 13)
    ax.set_yscale('log')
    ax.set_title('2D plot completness', fontsize = 13)
    fig.colorbar(c, ax=ax, label = 'Completness')
    plt.savefig(outpath + 'completeness_2D_plot' + '.png')
    plt.close()
