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
matching_folder = matching_folder + 'amico_cosmoDC2/'

##########select case
catalog1 = 'c1.fits'
catalog2 = 'c2.fits'
##########

#matching = 'redMaPPer cosmoDC2 small: $\Lambda>0$, $m_{halo}>10^{13}$'
#matching = 'WaZP cosmoDC2 small: NGALS>0, $m_{halo}>10^{13}$'
matching = 'AMICO cosmoDC2 small: $m_{halo}>10^{13}$'

outpath = "/pbs/home/t/tguillem/web/clusters/cluster_challenge/debug/amico_cosmoDC2/"
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
#zbins = np.linspace(0.2,1.5,14)
zbins = np.linspace(0.2,1.2,11)
mbins = [10**13,10**13.5,10**14,10**14.5,10**15]
#zbins = np.linspace(0, 1.6, 9)
#mbins = np.logspace(14, 15, 5)
fig = plt.figure()#figsize=(figx,figy))
info = r_cf.plot(c2, col1='z', col2='mass', bins1=zbins, bins2=mbins, matching_type='cross', legend_format=lambda x: f'10^{{{np.log10(x)}}}', lines_kwargs_list = [{'color':'black'}, {'color':'red'}, {'color':'blue'}, {'color':'purple'}])
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
mbins = np.linspace(13, 15, 17)
plt.figure()
info = r_cf.plot2D(c2, col1='z', col2='log_mass', bins1=zbins, bins2=mbins,
                   matching_type='cross', plt_kwargs={'cmap':'jet'})#{'cmap':'nipy_spectra'} #jet #gnuplot
plt.savefig(outpath+'recovery_plot2D.png', bbox_inches='tight') 

#purity
zbins = np.linspace(0.2,1.5,14)
mbins = [0,10,20,30,100]
#zbins = np.linspace(0, 1.6, 9)
#mbins = np.logspace(14, 15, 5)
fig = plt.figure()
info = r_cf.plot(c1, col1='z', col2='mass', bins1=zbins, bins2=mbins, matching_type='cross', lines_kwargs_list = [{'color':'black'}, {'color':'red'}, {'color':'blue'}, {'color':'purple'}])
info['ax'].set_xlabel('$z_{cl}$')
info['ax'].set_ylabel('Purity') 
info['ax'].set_ylim(0,1.2)
info['ax'].set_xlim(0.2,1.6) 
info['ax'].set_title(matching)
plt.savefig(outpath+'purity_plot.png', bbox_inches='tight')
plt.close(fig)

#recovery_plot2D
plt.figure()
info = r_cf.plot2D(c1, col1='z', col2='mass', bins1=zbins, bins2=mbins,
                   matching_type='cross', plt_kwargs={'cmap':'jet'})
plt.savefig(outpath+'purity_plot2D.png', bbox_inches='tight') 
#sys.exit()

###CHECK: ra/dec map
plt.figure()
plt.plot(c1.data['ra'],c1.data['dec'],'rx', color = 'red', label = 'AMICO')
plt.plot(c2.data['ra'],c2.data['dec'],'b.', color = 'blue', label = 'cosmoDC2')
plt.xlim([60, 72])
plt.ylim([-50, -30])
plt.xlabel("ra")
plt.ylabel("dec");
plt.title(matching)
plt.legend(bbox_to_anchor = (1, 1), loc = 'upper right', prop = {'size': 15})
plt.savefig(outpath+"clusters.png", bbox_inches='tight')
plt.close()
