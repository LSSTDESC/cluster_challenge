
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

#linestyle = ('-','--') if cats[0][0]==cats[1][0] else ('-','-')

output_matched_catalog(inpath + '.'.join(cats[0]) + '.fits', inpath + '.'.join(cats[1]) + '.fits', inpath + 'output_catalog_12.fits', cl[0], cl[1], matching_type='cross', overwrite=True)
c_merged_12 = ClCatalog.read(inpath+'output_catalog_12.fits', name='merged', full=True)

matching = 'AMICO cosmoDC2: $m_{200c}>10^{14}$, NGALS>0'
labels=['z incl.']
colors=['black','red','blue','purple', 'cyan', 'green', 'brown', 'orange', 'indigo']
zlabels = ['0-0.2', '0.2-0.5','0.5-0.8','0.8-1.0','1.0-1.2','1.2-1.5','1.5-1.8']

nbins_x = 9
snrbins = [3,3.5,4,4.5,5,6,7,8,9]
zbins = [0,0.2,0.5,0.8,1.0,1.2,1.5,1.8]
nbins_z = len(zbins)

#prepare tables

c_merged=c_merged_12[c_merged_12['cat1_log10_mass']<100]
#c_merged=c_merged[c_merged['cat2_z_cl']<1.15]
#c_merged=c_merged[c_merged['cat1_z_cl']<1.15]
c_merged_cut=c_merged[c_merged['cat1_log10_mass']>14]
                      
c_halos = cl[0][cl[0]['log10_mass']>13]
#c_halos = c_halos[c_halos['z_cl']<1.15]
c_halos_cut = c_halos[c_halos['log10_mass']>14]

c_clusters = cl[1][cl[1]['snr_cl']>0]
#c_clusters = c_clusters[c_clusters['z_cl']<1.15]

compl_snr = np.empty([nbins_z, nbins_x])

bin_z = np.empty([nbins_z-1])
for ix in range(nbins_z-1):
    bin_z[ix] = 0.5 * (zbins[ix] + zbins[ix+1])

for i in range(0,len(zbins)-1):
    zcut1 = zbins[i]
    zcut2 = zbins[i+1]
    c_halos_matched = c_merged_cut[(c_merged_cut['cat1_z_cl']>zcut1)*(c_merged_cut['cat1_z_cl']<zcut2)]
    for j in range(0,nbins_x):
        print('-----'+str(j))
        cut1 = snrbins[j]
        c_halos_matched = c_halos_matched[(c_halos_matched['cat2_snr_cl']>cut1)]
        n_halos_matched = len(c_halos_matched)
        n_halos = len(c_halos_cut[(c_halos_cut['z_cl']>zcut1)*(c_halos_cut['z_cl']<zcut2)])
        print(n_halos_matched)
        print(n_halos)
        if n_halos!=0:
            compl_snr[i][j] = round(n_halos_matched/n_halos,4)
        else : 
            compl_snr[i][j] = np.nan
        print(compl_snr[i][j])
    plt.scatter(snrbins, compl_snr[i], label=zlabels[i], color=colors[i], marker= ".", s=30)
    plt.plot(snrbins, compl_snr[i], color=colors[i])

plt.ylim(0, 1.2)
plt.xlim(2,10)
plt.xlabel('SNR')
plt.ylabel('Completeness')
plt.title(matching)
plt.legend()
plt.savefig(outpath+"completeness_snr_zbins_m>14.png", bbox_inches='tight')
plt.close()

#purity versus SNR
#nbins_x = 8
#snrbins = [3,4,5,6,7,8,9,10,11]

purity_snr = np.empty([nbins_z, nbins_x])
for i in range(0,len(zbins)-1):
    zcut1 = zbins[i]
    zcut2 = zbins[i+1]
    c_clusters_matched = c_merged_cut[(c_merged_cut['cat2_z_cl']>zcut1)*(c_merged_cut['cat2_z_cl']<zcut2)]
    for j in range(0,nbins_x):
        print('-----'+str(j))
        cut1 = snrbins[j]
        c_clusters_matched = c_clusters_matched[c_clusters_matched['cat2_snr_cl']>cut1]
        c_clusters_all = c_clusters[(c_clusters['snr_cl']>cut1) * (c_clusters['z_cl']>zcut1)*(c_clusters['z_cl']<zcut2)]
        n_clusters_matched = len(c_clusters_matched)
        n_clusters_all = len(c_clusters_all)
        print(n_clusters_matched)
        print(n_clusters_all)
        if n_clusters_all!=0:
            purity_snr[i][j] = round(n_clusters_matched/n_clusters_all,4)
        else : 
            purity_snr[i][j] = np.nan
        print(purity_snr[i][j])
    plt.scatter(snrbins, purity_snr[i], label=zlabels[i], color=colors[i], marker= ".", s=30)
    plt.plot(snrbins, purity_snr[i], color=colors[i])
plt.ylim(0, 1.2)
plt.xlim(2,10)
plt.xlabel('SNR')
plt.ylabel('Purity')
plt.title(matching)
plt.legend()
plt.savefig(outpath+"purity_snr_zbins_m>14.png", bbox_inches='tight')
plt.close()

#completeness versus purity
for i in range(0,len(zbins)-1):
    plt.scatter(compl_snr[i], purity_snr[i], label=zlabels[i], color=colors[i], marker= ".", s=30)
    plt.plot(compl_snr[i], purity_snr[i], color=colors[i])
    plt.ylim(0, 1)
    plt.xlim(0, 1)
    plt.xlabel('Completeness')
    plt.ylabel('Purity')
    plt.title(matching)
    #plt.legend()
    print(snrbins, compl_snr, purity_snr)
    #for j, label in enumerate(snrbins):
    #    plt.annotate(str(label), (compl_snr[i][j], purity_snr[i][j]))
plt.legend()
plt.savefig(outpath+"purity_versus_completeness_zbins_m>14.png", bbox_inches='tight')

plt.close()

print('DONE')
sys.exit()
