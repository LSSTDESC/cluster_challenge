

import sys, yaml
import numpy as np

from clevar.catalog import ClCatalog

import matplotlib.pyplot as plt
from src.plot_style import set_style, figsize
import src.saving_figures as sfigs
import src.fit_func as ff
from src.pltFunc_scaling import *


# read config files as online arguments 
config = sys.argv[1]

# open config files
with open(config) as fstream:
	cfg = yaml.safe_load(fstream)

cat1 = cfg['cats']['cat1']
cat2 = cfg['cats']['cat2']
mt_method = cfg['matching']['method']
mt_pref = cfg['matching']['pref']
if mt_method == 'member' :
	mt_params = [cfg['matching']['minimum_share_fraction']]
else :
	mt_params = [cfg['matching']['delta_z'], cfg['matching']['match_radius']]

inpath = sfigs.make_path(cat1, cat2, mt_method, mt_pref, mt_params, base=cfg['paths']['performance']['in'])
outpath = sfigs.make_path(cat1, cat2, mt_method, mt_pref, mt_params, base=cfg['paths']['performance']['out'], addon='scaling')



## Create a merged catalog for the cross-matched pairs if it does not already exist.
cltags = cfg['merged_keys']
mtcats = ClCatalog.read(f"{inpath}merged_catalog.fits", 'merged', tags=cltags)


## Is this comparison between cosmoDC2 and a cl. finder or between two cl. finders?
if np.any([c.split('.')[0] == 'cosmoDC2' for c in [cat1, cat2]]) :
	compared_to_truth = True
else :
	compared_to_truth = False



if compared_to_truth :
	## richness vs mass
	x = mtcats['cat1_log10_mass']
	y = mtcats['cat2_log10_mass']
	yerr = ff.log_errors(mtcats['cat2_mass'], mtcats['cat2_mass_err'])
	xlabel = cfg['latex']['logmass']
	ylabel = f"$\log {cfg['latex']['richness_2'][1:]}"
	title = f"{'.'.join(cat1.split('.')[:-1])} - {'.'.join(cat2.split('.')[:-1])}"
	saveas = 'richness_mass'
	plot_richness_mass(x, y, yerr, masscut=13.5, xlabel=xlabel, ylabel=ylabel, title=title,
		fit='linear', outpath=outpath, saveas=saveas+'_linear')
	plot_richness_mass(x, y, yerr, xlabel=xlabel, ylabel=ylabel, title=title,
		fit='lawnchair', outpath=outpath, saveas=saveas+'_lawnchair')
	
	
	## richness vs mass in redshift bins
	Nzbins = 3
	zbins = np.linspace(0, max(mtcats['z_halo']), Nzbins+1)
	
	saveas = 'richness_mass_zbinned'
	
	## with linear fitting
	plot_richness_mass_zbinned(x, y, mtcats['z_halo'], yerr, zbins=zbins, xlabel=xlabel, ylabel=ylabel,
		fit='linear', outpath=outpath, saveas=saveas+'_linear')
	plot_richness_mass_zbinned(x, y, mtcats['z_halo'], yerr, masscut=13.5, zbins=zbins, xlabel=xlabel, ylabel=ylabel,
		fit='linear', outpath=outpath, saveas=saveas+'_linear_masscut')
	
	## with lawnchair fitting
	plot_richness_mass_zbinned(x, y, mtcats['z_halo'], yerr, zbins=zbins, xlabel=xlabel, ylabel=ylabel,
		fit='lawnchair', outpath=outpath, saveas=saveas+'_lawnchair')
	plot_richness_mass_zbinned(x, y, mtcats['z_halo'], yerr, masscut=13.5, zbins=zbins, xlabel=xlabel, ylabel=ylabel,
		fit='lawnchair', outpath=outpath, saveas=saveas+'_lawnchair_masscut')
	
	
	## richness Prob. binned in redshift and mass
	bin_fits = plot_richnessLogProb_binned(y, mtcats['z_halo'], mtcats['log_mass'], cfg,
		outpath=outpath, saveas='richnessProb_binned')
	
	## evolution of binned values with mass and richness
	plot_richnessLogProb_binned_popt(bin_fits, bintype='zbin', outpath=outpath, saveas='richnessProb_fits_zbinned',)
	plot_richnessLogProb_binned_popt(bin_fits, bintype='Mbin', outpath=outpath, saveas='richnessProb_fits_Mbinned',)

else:
	x = np.log10(mtcats['mass'])
	y = np.log10(mtcats['richness'])
	yerr = ff.log_errors(mtcats['richness'], mtcats['richness_err'])
	xlabel = f"$\log {cfg['latex']['mass'][1:]}"
	ylabel = f"$\log {cfg['latex']['richness_2'][1:]}"
	title = f"{'.'.join(cat1.split('.')[:-1])} - {'.'.join(cat2.split('.')[:-1])}"
	saveas = 'richness_richness'
	plot_richness_richness(x, y, yerr, richcut=-np.inf, xlabel=xlabel, ylabel=ylabel, label=True, title=title,
		fit='linear', outpath=outpath, saveas=saveas)
	
	## richness vs richness in redshift bins
	Nzbins = 2
	zbins = np.linspace(0, max(mtcats['z_halo']), Nzbins+1)
	
	saveas = 'richness_richness_zbinned'
	
	## with linear fitting
	plot_richness_richness_zbinned(x, y, mtcats['z_halo'], yerr, richcut=-np.inf, zbins=zbins, title=title, xlabel=xlabel, ylabel=ylabel, label=True, fit='linear', outpath=outpath, saveas=saveas)
	
