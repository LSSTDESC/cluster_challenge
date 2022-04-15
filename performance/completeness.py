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
from clevar.match_metrics.recovery import ClCatalogFuncs as r_cf

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
#matching = 'redMaPPer cosmoDC2 small: $\Lambda>0$, $m_{halo}>10^{13}$'
matching = 'WaZP cosmoDC2 small: NGALS>0, $m_{halo}>10^{13}$'

outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/debug/wazp_cosmoDC2/"
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

#plot style
figx=10
figy=7

#recovery_plot
zbins = np.linspace(0.2,1.5,14)
mbins = [10**13,10**13.5,10**14,10**14.5,10**15]
#zbins = np.linspace(0, 1.6, 9)
#mbins = np.logspace(14, 15, 5)
fig = plt.figure(figsize=(figx,figy))
info = r_cf.plot(c2, col1='z', col2='mass', bins1=zbins, bins2=mbins, matching_type='cross', legend_format=lambda x: f'10^{{{np.log10(x)}}}', lines_kwargs_list = [{'color':'black'}, {'color':'red'}, {'color':'blue'}, {'color':'purple'}])
info['ax'].set_xlabel('$z_{halo}$')
info['ax'].set_ylabel('Completeness') 
info['ax'].set_ylim(0,1.2)
info['ax'].set_xlim(0.2,1.6) 
info['ax'].set_title(matching)
plt.savefig(outpath+'recovery_plot.png', bbox_inches='tight')
plt.close(fig)

#recovery_plot_panel
plt.figure()
info = r_cf.plot_panel(c2, col1='z', col2='mass', bins1=zbins, bins2=mbins,
                       matching_type='cross', label_format=lambda x: f'10^{{{np.log10(x)}}}')
plt.savefig(outpath+'recovery_plot_panel.png', bbox_inches='tight')
plt.close()

#recovery_plot2D
plt.figure()
info = r_cf.plot2D(c2, col1='z', col2='mass', bins1=zbins, bins2=mbins,
                   matching_type='cross', scale2='log')
plt.savefig(outpath+'recovery_plot2D.png', bbox_inches='tight') 
sys.exit()
