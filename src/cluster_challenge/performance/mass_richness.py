

import sys, yaml
import numpy as np
from astropy.table import Table

from clevar.catalog import ClCatalog
from clevar.match import get_matched_pairs, output_matched_catalog
from clevar.match_metrics import scaling

import matplotlib.pyplot as plt
from plot_style import set_style, figsize
import saving_figures as sfigs
import fit_func as ff
from plotting_functions import plot_richness_mass, plot_richness_mass_zbinned, plot_richnessLogProb_binned, plot_richnessLogProb_binned_popt
from scipy.optimize import curve_fit


# read config files as online arguments 
config = sys.argv[1]

# open config files
with open(config) as fstream:
	cfg = yaml.safe_load(fstream)

inpath = sfigs.make_path(
	cat1 = cfg['cats']['cat1'],
	cat2 = cfg['cats']['cat2'],
	mt_method = cfg['matching']['method'],
	mt_pref = cfg['matching']['pref'],
	mt_params = cfg['matching']['params'],
	base = cfg['path_base']['in']
)

outpath = sfigs.make_path(
	cat1 = cfg['cats']['cat1'],
	cat2 = cfg['cats']['cat2'],
	mt_method = cfg['matching']['method'],
	mt_pref = cfg['matching']['pref'],
	mt_params = cfg['matching']['params'],
	base = cfg['path_base']['out'],
	addon = 'testing/mass_richness'
)


## Select the catalogs.
cat1 = cfg['cats']['cat1']
cat2 = cfg['cats']['cat2']


## Create a merged catalog for the cross-matched pairs if it does not already exist.
cltags = cfg['merged_keys']
try :
	mtcats = ClCatalog.read(f"{inpath}{cat1}_{cat2}_merged.fits", 'merged', tags=cltags)
except :
	output_matched_catalog(f"{inpath}{cat1}.fits", f"{inpath}{cat2}.fits", f"{inpath}{cat1}_{cat2}_merged.fits",
		cats[0], cats[1], matching_type='cross', overwrite=True)
	mtcats = ClCatalog.read(f"{inpath}{cat1}_{cat2}_merged.fits", 'merged', tags=cltags)


## richness vs mass
x = mtcats['log_mass']
y = np.log10(mtcats['richness'])
yerr = ff.log_errors(mtcats['richness'], mtcats['richness_err'])
xlabel = cfg['latex']['logmass']
ylabel = f"$\log {cfg['latex']['richness_2'][1:]}"
title = f"{'.'.join(cat1.split('.')[:-1])} - {'.'.join(cat2.split('.')[:-1])}"
saveas = 'richness_mass'
plot_richness_mass(x, y, yerr, masscut=13.5, xlabel=xlabel, ylabel=ylabel, label=True, title=title, fit='linear', outpath=outpath, saveas=saveas+'_linear')
plot_richness_mass(x, y, yerr, xlabel=xlabel, ylabel=ylabel, label=True, title=title, fit='lawnchair', outpath=outpath, saveas=saveas+'_lawnchair')


## richness vs mass in redshift bins
Nzbins = 3
zbins = np.linspace(0, max(mtcats['z_halo']), Nzbins+1)

saveas = 'richness_mass_zbinned'

## with linear fitting
plot_richness_mass_zbinned(x, y, mtcats['z_halo'], yerr, masscut=0, zbins=zbins, xlabel=xlabel, ylabel=ylabel, label=True, fit='linear', outpath=outpath, saveas=saveas+'_linear')
plot_richness_mass_zbinned(x, y, mtcats['z_halo'], yerr, masscut=13.5, zbins=zbins, xlabel=xlabel, ylabel=ylabel, label=True, fit='linear', outpath=outpath, saveas=saveas+'_linear_masscut')

## with lawnchair fitting
plot_richness_mass_zbinned(x, y, mtcats['z_halo'], yerr, masscut=0, zbins=zbins, xlabel=xlabel, ylabel=ylabel, label=True, fit='lawnchair', outpath=outpath, saveas=saveas+'_lawnchair')
plot_richness_mass_zbinned(x, y, mtcats['z_halo'], yerr, masscut=13.5, zbins=zbins, xlabel=xlabel, ylabel=ylabel, label=True, fit='lawnchair', outpath=outpath, saveas=saveas+'_lawnchair_masscut')


## richness Prob. binned in redshift and mass
bin_fits = plot_richnessLogProb_binned(y, mtcats['z_halo'], mtcats['log_mass'], cfg, outpath=outpath, saveas='richnessProb_binned')

## evolution of binned values with mass and richness
plot_richnessLogProb_binned_popt(bin_fits, bintype='zbin', outpath=outpath, saveas='richnessProb_fits_zbinned',)
plot_richnessLogProb_binned_popt(bin_fits, bintype='Mbin', outpath=outpath, saveas='richnessProb_fits_Mbinned',)
