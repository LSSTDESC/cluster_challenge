
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
from plotting_functions import *

from tqdm import tqdm



## READ IN USER INPUT
cats = [sys.argv[1].split('_'), sys.argv[2].split('_')]
mt_method = sys.argv[3]
mt_pref = sys.argv[4]
if mt_method == 'proximity' :
	mt_params = [float(sys.argv[5]), float(sys.argv[6])]
elif mt_method == 'member' :
	mt_params = [float(sys.argv[5])]


## MAKE COSMODC2 CAT AS FIRST INDEX; OTHERWISE, SORT ALPHABETICALLY
## (i.e. SORTING ORDER: COSMODC2, AMIC0, REDMAPPER, WAZP)
if np.any([c[0] == 'cosmoDC2' for c in cats]) :
	first = np.argwhere(np.array([c[0] == 'cosmoDC2' for c in cats]))[0][0]
	last  = np.argwhere(np.array([c[0] != 'cosmoDC2' for c in cats]))[0][0]
else :
	(first, last) = np.argsort([c[0] for c in cats])



area = [440 if 'cosmoDC2' in c else 300 if 'DC2' in c else None for c in (sys.argv[1], sys.argv[2])]



## EXAMPLE INPUT: cosmoDC2_v0 wazp_cosmoDC2.fzb_v0 proximity more_massive 0.05 1
## RESULTING INPATH: /cosmoDC2_wazp.cosmoDC2.fzb/v0_v0/proximity_matching/deltaz_0.05_matchradius_1.0mpc_pref_more_massive/
inpath  = make_inpath(cats, mt_method, mt_pref, mt_params)
outpath = make_outpath(cats, mt_method, mt_pref, mt_params)



## READ IN MATCHED CATALOGS
cl = []
for c in cats :
	cl.append(Table(fits.open(inpath + '.'.join(c) + '.fits')[1].data))


cl[0].sort('mt_cross')
cl[1].sort('id_cl')

matches = [cl[0]['mt_cross'] != '', cl[1]['mt_cross'] != '']


## USE CONSISTENT COLORS IN PLOTS
colors = {'cosmoDC2':'black',
	'cosmoDC2.fzb':'gray',
	'wazp':'C0',
	'redmapper':'red'}

linestyle = ('-','--') if cats[0][0]==cats[1][0] else ('-','-')

### ||||||||||||||||||||||||||||||[ PLOT ON SKY ]||||||||||||||||||||||||||||||

## PLOT ON SKY MAP OF THE CLFINDER CLUSTERS
if True :
	for i in range(len(cl)) :
		plot_on_sky(cl[i]['ra_cl'], cl[i]['dec_cl'], alpha=0.2, title=f"{'.'.join(cats[i][:-1])}",
			label=('halos' if (cats[i][0]=='cosmoDC2') else 'clusters'),
			outpath=outpath, saveas=f"{'.'.join(cats[i])}_CLUSTER_ONSKY")

	plot_on_sky(cl[0]['ra_cl'][matches[0]], cl[0]['dec_cl'][matches[0]], alpha=0.2,
		title=f"{'.'.join(cats[first][:-1])}-{'.'.join(cats[last][:-1])}",
		label=('cl-halo matches' if ('cosmoDC2' in (cats[0][0],cats[1][0])) else 'cl-cl matches'),
		outpath=outpath, saveas=f"{'.'.join(cats[0])}_{'.'.join(cats[1])}_MATCHED_CLUSTER_ONSKY")

## UPDATE INDEX.HTML FILE
update_index_file_for_html_display(outpath+'onsky/', f"Onsky maps for {'.'.join(cats[0])} and {'.'.join(cats[1])}")


## |||||||||||||||||||||||||||||[ PLOT REDSHIFT ]|||||||||||||||||||||||||||||

## PLOT REDSHIFT HISTOGRAMS
if True :
	fill=False
	if cats[0][0] == cats[1][0] :
		fill = True
	
	if np.any([c[0] == 'cosmoDC2' for c in cats]) :
		annotate = 'clusters/halos'
	else :
		annotate = 'clusters'
	
	fig = plt.subplots(1,1, figsize=(3.5,3.5))
	fig[1].ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	bin_width = 0.05
	plot_redshift_hist(cl[0]['z_cl'], bin_width=bin_width, fill=fill, color=colors[cats[0][0]], label=f"{'.'.join(cats[0][:-1])}",
	        only_layer=False, fig=fig, last_layer=False)
	plot_redshift_hist(cl[1]['z_cl'], bin_width=bin_width, fill=False, color=colors[cats[1][0]], label=f"{'.'.join(cats[1][:-1])}",
		annotate=annotate, only_layer=False, outpath=outpath, fig=fig,
		saveas=f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_CLUSTER_REDSHIFT_HIST")
	
	fig = plt.subplots(1,1, figsize=(3.5,3.5))
	fig[1].ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])


	plot_redshift_hist(cl[0]['z_cl'][matches[0]], bin_width=bin_width, fill=fill, color=colors[cats[0][0]],
		label=f"{'.'.join(cats[0][:-1])}", only_layer=False, fig=fig, last_layer=False)
	plot_redshift_hist(cl[1]['z_cl'][matches[1]], bin_width=bin_width, fill=False, color=colors[cats[1][0]],
		label=f"{'.'.join(cats[1][:-1])}", title=f'Matching: {mt_method}', annotate=annotate,  only_layer=False, fig=fig,
		outpath=outpath, saveas=f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_MATCHED_CLUSTER_REDSHIFT_HIST")

## PLOT Z VS Z HIST2D
if True :
	if np.any([c[0] == 'cosmoDC2' for c in cats]) :		## IF COSMODC2 IS USED WE WANT TRUEZ TO BE ALONG HORIZONTAL
		xs = cl[first]['z_cl'][matches[first]]
		xlabel = "$z_{true}$"
		ys = cl[last]['z_cl'][matches[last]]
		ylabel = f"$z_{{{'.'.join(cats[last][:-1])}}}$"
	else :
		xs = cl[0]['z_cl'][matches[0]]
		xlabel = f"$z_{{{'.'.join(cats[0][:-1])}}}$"
		ys = cl[1]['z_cl'][matches[1]]
		ylabel = f"$z_{{{'.'.join(cats[1][:-1])}}}$"
	plot_redshift_zVSz(xs, ys, xlabel=xlabel, ylabel=ylabel, title=f'Matching: {mt_method}', diagonal=True,
		outpath=outpath, saveas=f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_MATCHED_CLUSTER_Z_VS_Z")


## PLOT BIAS AND STD OF REDSHIFT OBS VS TRUE
if True :
	if np.any([c[0] == 'cosmoDC2' for c in cats]) :		## IF COSMODC2 IS USED WE WANT TRUEZ TO BE ALONG HORIZONTAL
		xs = cl[first]['z_cl'][matches[first]]
		xlabel = "$z_{true}$"
		ys = cl[last]['z_cl'][matches[last]]
		ylabel = f"$z_{{cats[last]}}$"
		title=f"{'.'.join(cats[first][:-1])}-{'.'.join(cats[last][:-1])}"

		plot_redshift_std_and_mean(xs, ys, title=title,
			outpath=outpath, saveas=f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_MATCHED_CLUSTER_REDSHIFT_METRICS")

## UPDATE INDEX.HTML FILE
update_index_file_for_html_display(outpath+'redshifts/', f"Redshift distributions for {'.'.join(cats[0])} and {'.'.join(cats[1])}")


## |||||||||||||||||||||||||||||[ PLOT RICHNESS ]||||||||||||||||||||||||||||


## PLOT CLUSTER RICHNESS_1 VS HALO MASS OR CLUSTER RICHNESS_2
if True :
	if np.any([c[0] == 'cosmoDC2' for c in cats]) :
		xs = cl[first]['log10_mass'][matches[first]]
		xlabel = '$\log M_{200c}$'
		ys = np.log10(cl[last]['mass'][matches[last]])
		yerrs = 1/cl[last]['snr_cl'][matches[last]]
		ylabel = f"$\log \lambda_{{{'.'.join(cats[last][:-1])}}}$"
		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_MATCHED_HALO_MASS_VS_CLUSTER_RICHNESS"
		title=f"{'.'.join(cats[first][:-1])}-{'.'.join(cats[last][:-1])}"
		plot_richness_mass(xs, ys, yerrs, xlabel=xlabel, ylabel=ylabel, label=True, title=title, fit=True, outpath=outpath, saveas=saveas)
	else :
		xs = np.log10(cl[first]['mass'][matches[first]])
		xlabel = f"$\log \lambda_{{{'.'.join(cats[first][:-1])}}}$"
		ys = np.log10(cl[last]['mass'][matches[last]])
		ylabel = f"$\log \lambda_{{{'.'.join(cats[last][:-1])}}}$"
		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_MATCHED_CLUSTER_RICHNESS_VS_RICHNESS"
		title=None
		plot_richness_richness(xs, ys, xlabel=xlabel, ylabel=ylabel, label=True, title=title, fit=True, outpath=outpath, saveas=saveas)


## PLOT CLUSTER RICHNESS VS HALO MASS BINNED BY REDSHIFT
if True :
	def richness_mass_relation(x, c0, c1, c2, c3) :
		return 0.5*c0*((x-c2)/c1*erf((x-c2)/c1) + 1/np.sqrt(np.pi)*np.exp(-((x-c2)/c1)**2) + (x-c2)/c1) + c3
	
	def richness_richness_relation(x, c0, c1) :
		return c0 + c1*x

	zbin_width = 0.3
	zbins = np.linspace(0, 3, int((3-0) / zbin_width + 1))
	when_zbin = np.digitize(cl[first]['z_cl'][matches[first]], zbins) - 1
	Nzbins = max(when_zbin)

	if np.any([c[0] == 'cosmoDC2' for c in cats]) :
		params = [curve_fit(richness_mass_relation, cl[first]['log10_mass'][matches[first]][when_zbin==i],
			np.log10(cl[last]['mass'][matches[last]][when_zbin==i]), sigma=1/cl[last]['snr_cl'][matches[last]][when_zbin==i],
			bounds=((0.1,0.1,12,0),(5,5,15,2))) for i in range(Nzbins)]
		perrs = [np.sqrt(np.diag(ps[1])) for ps in params]
		
		fig, axs = plt.subplots(1, Nzbins, figsize=(3.5*Nzbins,3.5), sharey=True, sharex=True);
		Nhex = 30
		for i in range(Nzbins) :
			x4plot = np.linspace(13,16,100)

			slope = params[i][0][0]/params[i][0][1]
			slope_err = np.sqrt((perrs[i][0]/params[i][0][1])**2 + (params[i][0][0]*perrs[i][1]/params[i][0][1]**2)**2)
			
			axs[i].plot(x4plot, richness_mass_relation(x4plot, *params[i][0]), lw=1, c='r');
			axs[i].plot(x4plot, params[i][0][0]/params[i][0][1]*(x4plot - params[i][0][2]) + params[i][0][3], lw=1, c='r', ls='--',
				label='$\lambda \propto M_{{200c}}^{{{{{:.2f}}} \pm {{{:.2f}}}}}$'.format(slope, slope_err));
			axs[i].hexbin(cl[first]['log10_mass'][matches[first]][when_zbin==i],
				np.log10(cl[last]['mass'][matches[last]][when_zbin==i]),
				gridsize=(int(Nhex*hex_w),int(Nhex*hex_h)), norm=mpl.colors.LogNorm(), edgecolors='none');
			axs[i].set_ylim(0.25,2.5);
			axs[i].set_xlim(13,15.5);
			if i == 0 :
			    axs[i].set_ylabel(f"$\log \lambda_{{{'.'.join(cats[last][:-1])}}}$", loc='top');
			axs[i].set_xlabel('$\log M_{200c}$', loc='right');
			axs[i].annotate(f'${{{zbins[i]:.1f}}} < z_{{cl}} \leq {{{zbins[i+1]:.1f}}}$',
				xy=(axs[i].axis()[0], axs[i].axis()[3]), xycoords='data', va='top', ha='left',
				bbox=dict(boxstyle='round', fc='w', lw=1));
			axs[i].legend(loc='lower right', frameon=True);

		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_MATCHED_CLUSTER_RICHNESS_VS_HALO_MASS_BINNED_REDSHIFT"
		save_figure(fig, outpath+'richness/', saveas)

	else :
		params = [curve_fit(richness_richness_relation,
			np.log10(cl[first]['mass'][matches[first]][when_zbin==i]),
			np.log10(cl[last]['mass'][matches[last]][when_zbin==i]),
			bounds=((0,-2),(2,10))) for i in range(Nzbins)]
		perrs = [np.sqrt(np.diag(ps[1])) for ps in params]
		
		fig, axs = plt.subplots(1, Nzbins, figsize=(3.5*Nzbins,3.5), sharey=True, sharex=True);
		Nhex = 30
		for i in range(Nzbins) :
			x4plot = np.linspace(0,3,3)

			axs[i].plot(x4plot, richness_richness_relation(x4plot, *params[i][0]), lw=1, c='r',
				label='$\lambda_| = \lambda_{{-}}^{{\left({:.2f} \pm {:.2f}\\right)}} 10^{{\left({:.2f} \pm {:.2f}\\right)}} $'.format(
					params[i][0][1],perrs[i][1],params[i][0][0], perrs[i][0]))
			axs[i].hexbin(np.log10(cl[first]['mass'][matches[first]][zs_in_bin]),
				np.log10(cl[last]['mass'][matches[last]][zs_in_bin]),
				gridsize=(int(Nhex*hex_w),int(Nhex*hex_h)), norm=mpl.colors.LogNorm(), edgecolors='none')
			axs[i].set_ylim(0.25,3);
			axs[i].set_xlim(0.5,3);
			if i == 0 :
			    axs[i].set_ylabel(f"$\log \lambda_{{{'.'.join(cats[last][:-1])}}}$", loc='top');
			axs[i].set_xlabel(f"$\log \lambda_{{{'.'.join(cats[first][:-1])}}}$", loc='right');
			axs[i].annotate(f'${{{zbins[i]:.1f}}} < z_{{cl}} \leq {{{zbins[i+1]:.1f}}}$',
				xy=(axs[i].axis()[0], axs[i].axis()[3]), xycoords='data', va='top', ha='left',
				bbox=dict(boxstyle='round', fc='w', lw=1));
			axs[i].legend(loc='lower right', frameon=True);

		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_MATCHED_CLUSTER_RICHNESS_VS_RICHNESS_BINNED_REDSHIFT"
		save_figure(fig, outpath+'richness/', saveas)


## UPDATE INDEX.HTML FILE
update_index_file_for_html_display(outpath+'richness/', f"Richness plots for {'.'.join(cats[0])} and {'.'.join(cats[1])}")


## |||||||||||||||||||||||||||||[ PLOT DENSITY ]||||||||||||||||||||||||||||



## DENSITY PLOTS
if True :
	if not np.any([c[0] == 'cosmoDC2' for c in cats]) :
		density = [[np.sum(cl[i]['mass'] > x) / area[i] for x in np.logspace(0,3,300)] for i in range(2)]
	
		start_i = np.array([np.nonzero(np.diff(density[0]))[0][0] for i in range(2)])
		last_i  = np.array([np.nonzero(density[i])[0][-1] for i in range(2)])

		for i in range(2) :
			fig = plt.subplots(1,1, figsize=(3.5,3.5))
			if np.any([c[0] == 'cosmoDC2' for c in cats]) :
				xs = np.logspace(13,17,300)
				xlabel = '$M_{200c}$'
				ylabel = 'Density [halos/deg$^2$]'
			else :
				xs = np.logspace(0,3,300)
				xlabel = f"$\lambda_{{{'.'.join(cats[i][:-1])}}}$"
				ylabel = r'Density [clusters/deg$^2$]'
	
			fig[1].plot(xs[start_i[i]:last_i[i]], density[i][start_i[i]:last_i[i]], lw=1, c=colors[cats[i][0]])
			fig[1].set_xscale('log')
			fig[1].set_yscale('log')
			fig[1].set_xlabel(xlabel, loc='right')
			fig[1].set_ylabel(ylabel, loc='top')
			
			saveas = f"{'.'.join(cats[i])}_DENSITY"
			save_figure(fig[0], outpath+'density/', saveas)

		fig = plt.subplots(1,1, figsize=(3.5,3.5))

		xs = np.logspace(0,3,300)
		xlabel = r'$\lambda$'
		ylabel = r'Density [clusters/deg$^2$]'

		fig[1].plot(xs[start_i[first]:last_i[first]], density[first][start_i[first]:last_i[first]], lw=1,
			c=colors[cats[first][0]], ls=linestyle[0], label=f"{'.'.join(cats[first][:-1])}")
		fig[1].plot(xs[start_i[last]:last_i[last]], density[last][start_i[last]:last_i[last]], lw=1, c=colors[cats[last][0]], ls=linestyle[1],
			label=f"{'.'.join(cats[last][:-1])}")
		fig[1].set_xscale('log')
		fig[1].set_yscale('log')
		fig[1].set_xlabel(xlabel, loc='right')
		fig[1].set_ylabel(ylabel, loc='top')
		if cats[0][0] == cats[1][0] :
			handlelength = 1.5
			handletextpad = 0.5
			visible = True
		else :
			handlelength = 0
			handletextpad = 0
			visible = False
		leg = fig[1].legend(labelcolor='mec', frameon=False, fontsize='small', loc='upper right',
			handlelength=handlelength, handletextpad=handletextpad)
		for item in leg.legendHandles :
			item.set_visible(visible)
		
		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_DENSITY"
		save_figure(fig[0], outpath+'density/', saveas)

## UPDATE INDEX.HTML FILE
update_index_file_for_html_display(outpath+'density/', f"Density plots for {'.'.join(cats[0])} and {'.'.join(cats[1])}")


## |||||||||||||||||||||||||||||[ PLOT COMPLETENESS/PURITY ]||||||||||||||||||||||||||||


if True :
	zbins = np.linspace(0.2,1.1,12)
	mbins = np.linspace(13,15,5)
	richbins = np.linspace(0,3,5)

	c1 = ClCatalog.read_full(inpath + '.'.join(cats[first]) + '.fits')
	c2 = ClCatalog.read_full(inpath + '.'.join(cats[last]) + '.fits')

	if np.any([c[0] == 'cosmoDC2' for c in cats]) :
		fig = plt.figure(figsize=(3.5,3.5));
		info = r_cf.plot(c1, col1='z', col2='log_mass', bins1=zbins, bins2=mbins, matching_type='cross', legend_format=lambda x: f'{x}',
			lines_kwargs_list = [{'color':'black'}, {'color':'red'}, {'color':'blue'}, {'color':'purple'}]);
		info['ax'].set_xlabel('$z_{halo}$', loc='right');
		info['ax'].set_ylabel('Completeness', loc='top');
		info['ax'].set_ylim(0,1.05);
		info['ax'].set_xlim(0,1.5);
		# plt.legend(loc='lower right', frameon=False);
		leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title='$\log M_{200c}$', alignment='right', handlelength=0, handletextpad=0)
		for item in leg.legendHandles :
		    item.set_visible(False)
		title=f"{'.'.join(cats[first][:-1])}-{'.'.join(cats[last][:-1])}"
		fig.text(0,0, title, va='baseline', ha='left', size='small');
		
		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_COMPLETENESS_VS_REDSHIFT"
		save_figure(fig, outpath+'completeness_purity/', saveas)

	if np.any([c[0] == 'cosmoDC2' for c in cats]) :
		fig = plt.figure(figsize=(3.5,3.5));
		info = r_cf.plot(c2, col1='z', col2='mass', bins1=zbins, bins2=10**richbins, matching_type='cross',
			legend_format=lambda x: f'{np.log10(x)}',
			lines_kwargs_list = [{'color':'black'}, {'color':'red'}, {'color':'blue'}, {'color':'purple'}]);
		info['ax'].set_xlabel('$z_{cl}$', loc='right');
		info['ax'].set_ylabel('Purity', loc='top');
		info['ax'].set_ylim(0,1.05);
		info['ax'].set_xlim(0,1.5);
		# plt.legend(loc='lower right', frameon=False);
		leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title='$\log \lambda$', alignment='right',
			handlelength=0, handletextpad=0)
		for item in leg.legendHandles :
		    item.set_visible(False)
		title=f"{'.'.join(cats[first][:-1])}-{'.'.join(cats[last][:-1])}"
		fig.text(0,0, title, va='baseline', ha='left', size='small');
		
		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_PURITY_VS_REDSHIFT"
		save_figure(fig, outpath+'completeness_purity/', saveas)

	if np.any([c[0] == 'cosmoDC2' for c in cats]) :
		fig = plt.figure(figsize=(5,5))
		info = r_cf.plot_panel(c1, col1='z', col2='log_mass', bins1=zbins, bins2=mbins, matching_type='multi_self',
			label_format=lambda x: f'{x}')
		
		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_COMPLETENESS_VS_REDSHIFT_PANELS"
		save_figure(fig, outpath+'completeness_purity/', saveas)

	if np.any([c[0] == 'cosmoDC2' for c in cats]) :
		fig = plt.figure(figsize=(5,5))
		info = r_cf.plot_panel(c2, col1='z', col2='mass', bins1=zbins, bins2=10**richbins, matching_type='multi_self',
			label_format=lambda x: f'{np.log10(x)}')
		
		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_PURITY_VS_REDSHIFT_PANELS"
		save_figure(fig, outpath+'completeness_purity/', saveas)

	if np.any([c[0] == 'cosmoDC2' for c in cats]) :
		fig = plt.figure(figsize=(3.5,3.5))
		info = r_cf.plot2D(c1, col1='z', col2='log_mass', bins1=zbins, bins2=np.linspace(min(mbins), max(mbins), 30), matching_type='cross')
		plt.xlabel('$z_{halo}$', loc='right')
		plt.ylabel('$\log M_{200c}$', loc='top')
		title = f"{'.'.join(cats[first][:-1])} - {'.'.join(cats[last][:-1])}"
		plt.text(-0.1,12.7, title, va='baseline', ha='left', size='small')

		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_COMPLETENESS_2D"
		save_figure(fig, outpath+'completeness_purity/', saveas)

	if np.any([c[0] == 'cosmoDC2' for c in cats]) :
		fig = plt.figure(figsize=(3.5,3.5))
		info = recovery.plot(c1, 'cross', zbins[::2], np.logspace(min(mbins), max(mbins), 30), shape='line', transpose=True)
		plt.xlabel('$M_{200c}$', loc='right');
		plt.ylabel('completeness', loc='top');
		leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title='$z_{halo}$', alignment='right', handlelength=0, handletextpad=0)
		for item in leg.legendHandles :
		    item.set_visible(False)
		title=f"{'.'.join(cats[first][:-1])}-{'.'.join(cats[last][:-1])}"
		fig.text(0,0, title, va='baseline', ha='left', size='small');
		
		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_COMPLETENESS_VS_M200C"
		save_figure(fig, outpath+'completeness_purity/', saveas)

	if np.any([c[0] == 'cosmoDC2' for c in cats]) :
		fig = plt.figure(figsize=(3.5,3.5))
		info = recovery.plot(c2, 'cross', zbins[::2], np.logspace(min(richbins), max(richbins), 30), shape='line', transpose=True)
		plt.xlabel('$\lambda$', loc='right');
		plt.ylabel('purity', loc='top');
		leg = plt.legend(labelcolor='mec', frameon=True, fontsize='small', title='$z_{halo}$', alignment='right', handlelength=0, handletextpad=0)
		for item in leg.legendHandles :
		    item.set_visible(False)
		title=f"{'.'.join(cats[first][:-1])}-{'.'.join(cats[last][:-1])}"
		fig.text(0,0, title, va='baseline', ha='left', size='small');
		
		saveas = f"{'.'.join(cats[first])}_{'.'.join(cats[last])}_PURITY_VS_RICHNESS"
		save_figure(fig, outpath+'completeness_purity/', saveas)

## UPDATE INDEX.HTML FILE
update_index_file_for_html_display(outpath+'completeness_purity/', f"Completeness and Purity plots for {'.'.join(cats[0])} and {'.'.join(cats[1])}")


