

import sys, yaml
import numpy as np
from astropy.table import Table

from clevar.catalog import ClCatalog
from clevar.match import get_matched_pairs, output_matched_catalog
from clevar.match_metrics import scaling, distances
from clevar.match_metrics.recovery import ClCatalogFuncs as r_cf
from clevar.cosmology import AstroPyCosmology

import matplotlib.pyplot as plt
from plot_style import set_style, figsize
import saving_figures as sfigs
#import fit_func as ff
from plotting_functions import make_bins
#from scipy.optimize import curve_fit


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
	base = cfg['paths']['performance']['in']
)

outpath = sfigs.make_path(
	cat1 = cfg['cats']['cat1'],
	cat2 = cfg['cats']['cat2'],
	mt_method = cfg['matching']['method'],
	mt_pref = cfg['matching']['pref'],
	mt_params = cfg['matching']['params'],
	base = cfg['paths']['performance']['out'],
)


## Select the catalogs.
cat1 = cfg['cats']['cat1']
cat2 = cfg['cats']['cat2']

cats = {cat1: ClCatalog.read_full(f"{inpath}{cat1}.fits"),
	cat2: ClCatalog.read_full(f"{inpath}{cat2}.fits")}


## create the cosmology object (matching cosmoDC2 simulation)
H0 = cfg['cosmology']['h'] * 100
Omega_dm0 = cfg['cosmology']['Omega_dm0']
Omega_b0  = cfg['cosmology']['Omega_b0']
Omega_k0  = cfg['cosmology']['Omega_k0']
cosmo = AstroPyCosmology(H0=H0, Omega_b0=Omega_b0, Omega_dm0=Omega_dm0, Omega_k0=Omega_k0)


deg_bins = make_bins('radius', cfg, units='degrees', grain='both')
Mpc_bins = make_bins('radius', cfg, units='Mpc', grain='both')
mass_bins   = make_bins('mass', cfg, grain='both')
z_bins	    = make_bins('redshift', cfg, grain='both')
deltaz_bins = make_bins('delta_redshift', cfg, log=True, grain='both')
rich_bins   = make_bins('richness', cfg, grain='both')


## completeness vs redshift
for mt_type in ['cross', 'multi_self'] :
	fig = plt.figure(figsize=figsize())
	info = r_cf.plot(
		cats[cat1],
		col1=cfg['cat_keys']['cat1']['z'],
		col2=cfg['cat_keys']['cat1']['logmass'],
		bins1=z_bins['fine'],
		bins2=mass_bins['coarse'],
		matching_type=mt_type,
		legend_format=lambda x: f'{x:.1f}',
		lines_kwargs_list = [{'color':'black'}, {'color':'red'}, {'color':'blue'}, {'color':'purple'}]);
	info['ax'].set_xlabel(cfg['latex']['redshift_1'], loc='right')
	info['ax'].set_ylabel('Completeness', loc='top')
	info['ax'].set_ylim(0,1.05)
	info['ax'].set_xlim(min(z_bins['coarse']), max(z_bins['coarse']))
	leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['logmass'],
		alignment='right', handlelength=0, handletextpad=0)
	for item in leg.legendHandles :
	    item.set_visible(False)
	
	fig.text(0,0, f'{mt_type.split("_")[0]} matching', va='baseline', ha='left');

saveas = 'completeness_redshift'
sfigs.save_figure(fig, outpath, saveas=saveas)


## purity vs redshift
for mt_type in ['cross', 'multi_self'] :
	fig = plt.figure(figsize=figsize())
	info = r_cf.plot(
		cats[cat2],
		col1=cfg['cat_keys']['cat2']['z'],
		col2=cfg['cat_keys']['cat2']['richness'],
		bins1=z_bins['fine'],
		bins2=rich_bins['coarse'],
		matching_type=mt_type,
		legend_format=lambda x: f'{x}',
		lines_kwargs_list = [{'color':'black'}, {'color':'red'}, {'color':'blue'}, {'color':'purple'}])
	info['ax'].set_xlabel(cfg['latex']['redshift_1'], loc='right')
	info['ax'].set_ylabel('Purity', loc='top')
	info['ax'].set_ylim(0,1.05)
	info['ax'].set_xlim(min(z_bins['coarse']), max(z_bins['coarse']))
	leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small',
		title=cfg['latex']['richness_2'], alignment='right', handlelength=0, handletextpad=0)
	for item in leg.legendHandles :
	    item.set_visible(False)
	
	fig.text(0,0, f'{mt_type.split("_")[0]} matching', va='baseline', ha='left')

saveas = 'purity_redshift'
sfigs.save_figure(fig, outpath, saveas=saveas)



## Match count vs. angular separation
## binned by mass
info = distances.central_position(
	cats[cat1], cats[cat2], 'cross', radial_bins=deg_bins['fine'], radial_bin_units='degrees',
	quantity_bins=cfg['cat_keys']['cat1']['mass'], bins=10**mass_bins['coarse'], log_quantity=True)

bin_positions = (info['data']['distance_bins'][1:] + info['data']['distance_bins'][:-1])/2

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(mass_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]
	label = f"$[10^{{{mass_bins['coarse'][i]:.1f}}}10^{{{mass_bins['coarse'][i+1]:.1f}}}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['mass'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)
axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlabel('$\Delta\\theta$ [deg]', loc='right')
axs.set_ylabel('Number of matches', loc='top')

saveas = 'matches_degrees_massBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



bin_widths = (info['data']['distance_bins'][1:] - info['data']['distance_bins'][:-1])

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(mass_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i] / bin_widths
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=f"$[10^{{{mass_bins['coarse'][i]:.1f}}}:10^{{{mass_bins['coarse'][i+1]:.1f}}}]$")
axs.set_ylabel('No. of matches / deg', loc='top')
axs.set_xlabel('$\Delta\\theta$ [deg]', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['mass'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerDeg_degrees_massBins'
sfigs.save_figure(fig, outpath, saveas=saveas)


## binned by redshift
info = distances.central_position(
	cats[cat1], cats[cat2], 'cross', radial_bins=deg_bins['fine'], radial_bin_units='degrees',
	quantity_bins=cfg['cat_keys']['cat1']['z'], bins=z_bins['coarse'], log_quantity=False)

bin_positions = (info['data']['distance_bins'][1:] + info['data']['distance_bins'][:-1])/2

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(z_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]
	label = f"$[{z_bins['coarse'][i]:.1f}:{z_bins['coarse'][i+1]:.1f}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['redshift_1'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)
axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlabel('$\Delta\\theta$ [deg]', loc='right')
axs.set_ylabel('Number of matches', loc='top')

saveas = 'matches_degrees_redshiftBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(z_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]/bin_widths
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=f"$[{z_bins['coarse'][i]:.1f}:{z_bins['coarse'][i+1]:.1f}]$")
axs.set_ylabel('No. of matches / deg', loc='top')
axs.set_xlabel('$\Delta\\theta$ [deg]', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['redshift_1'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerDeg_degrees_redshiftBins'
sfigs.save_figure(fig, outpath, saveas=saveas)


## binned by richness
info = distances.central_position(
	cats[cat2], cats[cat1], 'cross', radial_bins=deg_bins['fine'], radial_bin_units='degrees',
	quantity_bins=cfg['cat_keys']['cat2']['richness'], bins=rich_bins['coarse'], log_quantity=False)

bin_positions = (info['data']['distance_bins'][1:] + info['data']['distance_bins'][:-1])/2

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(rich_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]
	label = f"$[{rich_bins['coarse'][i]:.1f}:{rich_bins['coarse'][i+1]:.1f}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['richness_2'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)
axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlabel('$\Delta\\theta$ [deg]', loc='right')
axs.set_ylabel('Number of matches', loc='top')

saveas = 'matches_degrees_richBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(rich_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]/bin_widths
	label = f"$[{rich_bins['coarse'][i]:.1f}:{rich_bins['coarse'][i+1]:.1f}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
axs.set_ylabel('No. of matches / deg', loc='top')
axs.set_xlabel('$\Delta\\theta$ [deg]', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['richness_2'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerDeg_degrees_richBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



## Match count vs. radial separation [Mpc]

## binned by mass
info = distances.central_position(
	cats[cat1], cats[cat2], 'cross', radial_bins=Mpc_bins['fine'], radial_bin_units='Mpc', cosmo=cosmo,
	quantity_bins=cfg['cat_keys']['cat1']['mass'], bins=10**mass_bins['coarse'], log_quantity=True)

bin_positions = (info['data']['distance_bins'][1:] + info['data']['distance_bins'][:-1])/2

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(mass_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]
	label = f"$[10^{{{mass_bins['coarse'][i]:.1f}}}:10^{{{mass_bins['coarse'][i+1]:.1f}}}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['mass'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)
axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlabel('$\Delta R$ [Mpc]', loc='right')
axs.set_ylabel('Number of matches', loc='top')

saveas = 'matches_Mpc_massBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



bin_widths = (info['data']['distance_bins'][1:] - info['data']['distance_bins'][:-1])

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(mass_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i] / bin_widths
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=f"$[10^{{{mass_bins['coarse'][i]:.1f}}}:10^{{{mass_bins['coarse'][i+1]:.1f}}}]$")
axs.set_ylabel('No. of matches / Mpc', loc='top')
axs.set_xlabel('$\Delta R$ [Mpc]', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['mass'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerMpc_Mpc_massBins'
sfigs.save_figure(fig, outpath, saveas=saveas)


bin_areas = np.pi*(info['data']['distance_bins'][1:]**2 - info['data']['distance_bins'][:-1]**2)

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(mass_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i] / bin_areas
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=f"$[10^{{{mass_bins['coarse'][i]:.1f}}}:10^{{{mass_bins['coarse'][i+1]:.1f}}}]$")
axs.set_ylabel('No. of matches / Area [Mpc$^2$]', loc='top')
axs.set_xlabel('$\Delta R$ [Mpc]', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['mass'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerMpc2_Mpc_massBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



## binned by redshift
info = distances.central_position(
	cats[cat1], cats[cat2], 'cross', radial_bins=Mpc_bins['fine'], radial_bin_units='Mpc', cosmo=cosmo,
	quantity_bins=cfg['cat_keys']['cat1']['z'], bins=z_bins['coarse'], log_quantity=False)

bin_positions = (info['data']['distance_bins'][1:] + info['data']['distance_bins'][:-1])/2

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(z_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]
	label = f"$[{z_bins['coarse'][i]:.1f}:{z_bins['coarse'][i+1]:.1f}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['redshift_1'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)
axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlabel('$\Delta R$ [Mpc]', loc='right')
axs.set_ylabel('Number of matches', loc='top')

saveas = 'matches_Mpc_redshiftBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(z_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]/bin_widths
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=f"$[{z_bins['coarse'][i]:.1f}:{z_bins['coarse'][i+1]:.1f}]$")
axs.set_ylabel('No. of matches / Mpc', loc='top')
axs.set_xlabel('$\Delta R$ [Mpc]', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['redshift_1'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerMpc_Mpc_redshiftBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(z_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]/bin_areas
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=f"$[{z_bins['coarse'][i]:.1f}:{z_bins['coarse'][i+1]:.1f}]$")
axs.set_ylabel('No. of matches / Area [Mpc$^2$]', loc='top')
axs.set_xlabel('$\Delta R$ [Mpc]', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['redshift_1'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerMpc2_Mpc_redshiftBins'
sfigs.save_figure(fig, outpath, saveas=saveas)


## binned by richness
info = distances.central_position(
	cats[cat2], cats[cat1], 'cross', radial_bins=Mpc_bins['fine'], radial_bin_units='Mpc', cosmo=cosmo,
	quantity_bins=cfg['cat_keys']['cat2']['richness'], bins=rich_bins['coarse'], log_quantity=False)

bin_positions = (info['data']['distance_bins'][1:] + info['data']['distance_bins'][:-1])/2

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(rich_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]
	label = f"$[{rich_bins['coarse'][i]:.1f}:{rich_bins['coarse'][i+1]:.1f}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['richness_2'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)
axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlabel('$\Delta R$ [Mpc]', loc='right')
axs.set_ylabel('Number of matches', loc='top')

saveas = 'matches_Mpc_richBins'
sfigs.save_figure(fig, outpath, saveas=saveas)


fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(rich_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i] / bin_widths
	label = f"$[{rich_bins['coarse'][i]:.1f}:{rich_bins['coarse'][i+1]:.1f}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
axs.set_ylabel('No. of matches / Mpc', loc='top')
axs.set_xlabel('$\Delta R$ [Mpc]', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['richness_2'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerMpc_Mpc_richBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(rich_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i] / bin_areas
	label = f"$[{rich_bins['coarse'][i]:.1f}:{rich_bins['coarse'][i+1]:.1f}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
axs.set_ylabel('No. of matches / Area [Mpc$^2$]', loc='top')
axs.set_xlabel('$\Delta R$ [Mpc]', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['richness_2'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerMpc2_Mpc_richBins'
sfigs.save_figure(fig, outpath, saveas=saveas)


## Match count vs. redshift separation

## binned by mass
info = distances.redshift(
	cats[cat1], cats[cat2], 'cross', redshift_bins=deltaz_bins['fine'],
	quantity_bins=cfg['cat_keys']['cat1']['mass'], bins=10**mass_bins['coarse'], log_quantity=True, normalize='cat1')

bin_positions = (info['data']['distance_bins'][1:] + info['data']['distance_bins'][:-1])/2

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(mass_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]
	label = f"$[10^{{{mass_bins['coarse'][i]:.1f}}}:10^{{{mass_bins['coarse'][i+1]:.1f}}}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['mass'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)
axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlabel('$\Delta z / (1+z)$', loc='right')
axs.set_ylabel('Number of matches', loc='top')

saveas = 'matches_redshift_massBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



bin_widths = (info['data']['distance_bins'][1:] - info['data']['distance_bins'][:-1]) / (1 + bin_positions)

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(mass_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i] / bin_widths
	label=f"$[10^{{{mass_bins['coarse'][i]:.1f}}}:10^{{{mass_bins['coarse'][i+1]:.1f}}}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
axs.set_ylabel('No. of matches$\\times (1+z) / \Delta z$', loc='top')
axs.set_xlabel('$\Delta z / (1+z)$', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['mass'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerDeltaz_redshift_massBins'
sfigs.save_figure(fig, outpath, saveas=saveas)


## binned by redshift
info = distances.central_position(
	cats[cat1], cats[cat2], 'cross', radial_bins=deltaz_bins['fine'],
	quantity_bins=cfg['cat_keys']['cat1']['z'], bins=z_bins['coarse'], log_quantity=False)

bin_positions = (info['data']['distance_bins'][1:] + info['data']['distance_bins'][:-1])/2

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(z_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]
	label = f"$[{z_bins['coarse'][i]:.1f}:{z_bins['coarse'][i+1]:.1f}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['redshift_1'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)
axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlabel('$\Delta z / (1+z)$', loc='right')
axs.set_ylabel('Number of matches', loc='top')

saveas = 'matches_redshift_redshiftBins'
sfigs.save_figure(fig, outpath, saveas=saveas)



fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(z_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]/bin_widths
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=f"$[{z_bins['coarse'][i]:.1f}:{z_bins['coarse'][i+1]:.1f}]$")
axs.set_ylabel('No. of matches$\\times (1+z) / \Delta z$', loc='top')
axs.set_xlabel('$\Delta z / (1+z)$', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['redshift_1'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerDeltaz_redshift_redshiftBins'
sfigs.save_figure(fig, outpath, saveas=saveas)


## binned by richness
info = distances.redshift(
	cats[cat2], cats[cat1], 'cross', redshift_bins=deltaz_bins['fine'],
	quantity_bins=cfg['cat_keys']['cat2']['richness'], bins=rich_bins['coarse'], log_quantity=False, normalize='cat2')

bin_positions = (info['data']['distance_bins'][1:] + info['data']['distance_bins'][:-1])/2

fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(rich_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i]
	label = f"$[{rich_bins['coarse'][i]:.1f}:{rich_bins['coarse'][i+1]:.1f}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['richness_2'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)
axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlabel('$\Delta z / (1+z)$', loc='right')
axs.set_ylabel('Number of matches', loc='top')

saveas = 'matches_redshift_richBins'
sfigs.save_figure(fig, outpath, saveas=saveas)


fig, axs = plt.subplots(1,1, figsize=figsize())
for i in range(len(rich_bins['coarse'])-1) :
	bin_heights = info['data']['hist'][:,i] / bin_widths
	label = f"$[{rich_bins['coarse'][i]:.1f}:{rich_bins['coarse'][i+1]:.1f}]$"
	axs.step(bin_positions, bin_heights, where='mid', lw=1, label=label)
axs.set_ylabel('No. of matches$\\times (1+z) / \Delta z$', loc='top')
axs.set_xlabel('$\Delta z / (1+z)$', loc='right')

leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title=cfg['latex']['richness_2'], alignment='right', handlelength=0, handletextpad=0)
for item in leg.legendHandles :
    item.set_visible(False)

axs.grid(which='major', lw=0.5)
axs.grid(which='minor', lw=0.2)
axs.set_xscale('log')
axs.set_yscale('log')

saveas = 'matchesPerDeltaz_redshift_richBins'
sfigs.save_figure(fig, outpath, saveas=saveas)
