
import numpy as np
import glob
import os
import pickle
import struct

from scipy.stats import binned_statistic, trim_mean
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.optimize import minimize, curve_fit
from scipy.special import erf

from src.plot_style import set_style, figsize
import src.saving_figures as sfigs
from src.fit_func import linear, lawnchair, gauss


## Initialize style for plots.
set_style()


def make_inpath(cats, mt_method, mt_pref, mt_params, base='/sps/lsst/groups/clusters/cluster_comparison_project/after_matching/') :

	(first, last) = np.argsort([c[0] for c in cats])

	inpath  = base
	inpath += f"{'.'.join(cats[first][:-1])}_{'.'.join(cats[last][:-1])}/"          ## EXAMPLE INPUT: cosmoDC2_v1 wazp_cosmoDC2.fzb_v0 proximity
	inpath += f"{cats[first][-1]}_{cats[last][-1]}/"                                ## RESULTING INPATH: /cosmoDC2_wazp.cosmoDC2.fzb/v1_v0/
	
	if mt_method == 'proximity' :
	        inpath += f"proximity_matching/deltaz_{mt_params[0]}_matchradius_{mt_params[1]}mpc_pref_{mt_pref}/"
	else :
	        inpath += f"member_matching/fshare_{mt_params[0]}_pref_{mt_pref}/"
	
	return inpath


def make_outpath(cats, mt_method, mt_pref, mt_params, base='/sps/lsst/users/rsolomon/web/desc/cluster_comparison_project/') :

	(first, last) = np.argsort([c[0] for c in cats])

	outpath  = base
	outpath += f"{'.'.join(cats[first][:-1])}_{'.'.join(cats[last][:-1])}/"
	outpath += f"{cats[first][-1]}_{cats[last][-1]}/"
	
	if mt_method == 'proximity' :
	        outpath += f"proximity_matching/deltaz_{mt_params[0]}_matchradius_{mt_params[1]}mpc_pref_{mt_pref}/"
	else :
	        outpath += f"member_matching/fshare_{mt_params[0]}_pref_{mt_pref}/"

	return outpath


def make_bins(param, cfg, log=False, units=None, grain='both') :
	'''
	Creates log or linear  bins in the parameter of choice from the config file.
	***
	param: str
		name of parameter to bin (should match naming in config file)
	cfg: dict
		dict containing the range and No. for coarse and finely grained bins
	log: bool (default: False)
		use linear (False) or logarithmic (True) bin sizes
	grain: str (default: both)
		 coarse or fine binning

	Returns: dict of arrays, array
		If grain='both' then a dict containing a coarse and fine grain binning is returned.
		Else, an array of the bins is returned.
	'''

	if param == 'radius' :
		info = cfg['match_metrics']['bins'][param][units]
	else :
		info = cfg['match_metrics']['bins'][param]

	left  = info['range'][0]
	right = info['range'][1]
	
	def space(l, r, N) :
		if not log :
			return np.linspace(l, r, N)
		else :
			return np.logspace(l, r, N)

	if grain == 'both' :
		Ncoarse = info['N']['coarse']
		Nfine   = info['N']['fine']
		
		coarse_bins = space(left, right, Ncoarse)
		fine_bins = space(left, right, Nfine)

		return {'coarse':coarse_bins, 'fine':fine_bins}
	else :
		N = info['N'][grain]
		
		bins = space(left, right, N)

		return bins



def plot_hist(x, bins, xlabel, ylabel='Count', fill=False, color='C0',
	title=None, label=None, annotate=None, xscale='linear', yscale='linear',
	outpath=None, saveas=None, show=False,
	only_layer=True, fig=None, last_layer=True) :
	
	if fig is None :        ## IN CASE OF MULTIPLE DATASETS ON SAME PLOT
		fig, axs = plt.subplots(1,1, figsize=figsize())
	else :
		fig, axs = fig

	axs.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	if fill == True :
	        alpha = 0.3
	else :
	        alpha = 1

	h = axs.hist(x, bins=bins, histtype='step', fill=fill, color=color, alpha=alpha, label=label)
	
	if last_layer and label != None :
		#axs.set_ylim(top=1.3*np.max(axs.get_ylim()))
		leg = axs.legend(labelcolor='mec', frameon=True, fontsize='small', loc='lower left', handlelength=0, handletextpad=0)
		for item in leg.legendHandles :
			item.set_visible(False)
	if last_layer and annotate != None :
		axs.annotate(annotate, xy=(axs.axis()[1], axs.axis()[3]), xycoords='data', va='top', ha='right',
			bbox=dict(boxstyle='round', fc='w', lw=1), zorder=np.inf)
	if last_layer :
		axs.set_xscale(xscale)
		axs.set_yscale(yscale)
		axs.set_xlabel(xlabel, loc='right')
		axs.set_ylabel(ylabel, loc='top')
		if title != None :
			fig.text(0,0, title, va='baseline', ha='left', size='x-small')
		if (outpath != None) & (saveas != None) :
			sfigs.save_figure(fig, outpath, saveas)
		if show :
			plt.show()



def plot_on_sky(ra, dec, Nhex, title=None, label=None, outpath=None, saveas=None, show=False,
	only_layer=True, axs=None, last_layer=True) :
	
	if only_layer :         ## IN CASE OF MULTIPLE DATASETS, PLACE fig, axs = ... BEFORE FUNCTION CALL
	        fig, axs = plt.subplots(1,1, figsize=figsize())
	
	axs.hexbin(ra, dec, gridsize=Nhex, norm=mpl.colors.LogNorm());
	
	if title != None :
		fig.text(0,0, title, va='baseline', ha='left', size='x-small')
	if last_layer and label != None :
		axs.annotate(label, xy=(axs.axis()[1], axs.axis()[3]), xycoords='data', va='top', ha='right', bbox=dict(boxstyle='round', fc='w', lw=1))
	if last_layer :
		axs.set_xlabel('RA', loc='right')
		axs.set_ylabel('DEC', loc='top')
		if saveas != None :
			sfigs.save_figure(fig, outpath, saveas)
		if show :
			plt.show()




def plot_redshift_redshift(z1, z2, Nhex=100, xlabel=None, ylabel=None, diagonal=False, title=None, outpath=None, saveas=None, show=False,) :
	fig, axs = plt.subplots(1,1, figsize=figsize())
	axs.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	if diagonal :
		axs.plot([min(z1),max(z1)],[min(z1),max(z1)], lw=0.5, c='r', ls='--')

	axs.hexbin(z1, z2, gridsize=Nhex, norm=mpl.colors.LogNorm());
	axs.grid(which='major', axis='both', ls='-', lw=0.5, color='grey')
	
	if title != None :
		fig.text(0,0, title, va='baseline', ha='left', size='x-small')
	        
	axs.set_xlabel(xlabel, loc='right')
	axs.set_ylabel(ylabel, loc='top')
	axs.set_xlim(0,1.5)
	axs.set_ylim(0,1.5)
	if saveas != None :
		sfigs.save_figure(fig, outpath, saveas, png=True, pdf=False, pkl=False)
	if show :
		plt.show()


def plot_zscatter_redshift(z1, z2, cfg, Nhex=100,
	label=None, title=None, outpath=None, saveas=None, show=False,
	only_layer=True, fig=None, last_layer=True) :

	fig = plt.subplots(2,1, figsize=figsize(2,1), gridspec_kw={'hspace':0.4}, sharex=True)
	[ax.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2]) for ax in fig[1]]
	[ax.grid(which='major', axis='both', ls='-', lw=0.5, color='grey') for ax in fig[1]]

	fig[1][0].hexbin(z1, (z2-z1), gridsize=Nhex, norm=mpl.colors.LogNorm())
	fig[1][1].hexbin(z2, (z2-z1), gridsize=Nhex, norm=mpl.colors.LogNorm())

	if title != None :
		fig[0].text(0,0, title, va='baseline', ha='left', size='x-small')
	        
	fig[1][0].set_xlabel(cfg['latex']['redshift_1'], loc='right')
	fig[1][1].set_xlabel(cfg['latex']['redshift_2'], loc='right')
	[ax.set_ylabel('$z_{cl} - z_{halo}$', loc='top') for ax in fig[1]]
	
	if saveas != None :
		sfigs.save_figure(fig[0], outpath, saveas, png=True, pdf=False, pkl=False)
	if show :
		plt.show()
	
	


def plot_redshift_std_and_mean(z1, z2, x_true=True, Nbins=50,
	label=None, title=None, outpath=None, saveas=None, show=False,
	only_layer=True, fig=None, last_layer=True) :

	zbins = np.linspace(min(z1), max(z1), Nbins)
	if x_true :
		vals = (z1-z2)/(1+z1)
		mean = binned_statistic(z1[vals<0.15], vals[vals<0.15], bins=Nbins, statistic='mean')[0]
		std  = binned_statistic(z1[vals<0.15], vals[vals<0.15], bins=Nbins, statistic='std')[0]
		xlabel = '$z_{halo}$'
	else :
		vals = (z1-z2)/(1+z2)
		mean = binned_statistic(z2[vals<0.15], vals[vals<0.15], bins=Nbins, statistic='mean')[0]
		std  = binned_statistic(z2[vals<0.15], vals[vals<0.15], bins=Nbins, statistic='std')[0]
		xlabel = '$z_{cl}$'

	if only_layer :
		fig = plt.subplots(2,1, figsize=figsize(2,1), gridspec_kw={'hspace':0.4}, sharex=True)
		[ax.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2]) for ax in fig[1]]
	
	fig[1][0].grid(which='major', axis='both', ls='-', lw=0.5, color='grey')
	fig[1][0].plot(zbins+0.5*np.mean(np.diff(zbins)), mean, lw=1)
	fig[1][0].axhline(trim_mean(mean, proportiontocut=0.1), ls='--', lw=1, c='k', label=f'$b_{{avg}} = ${trim_mean(mean, proportiontocut=0.1):.4f}')
	
	fig[1][1].grid(which='major', axis='both', ls='-', lw=0.5, color='grey')
	fig[1][1].plot(zbins+0.5*np.mean(np.diff(zbins)), std, lw=1)
	fig[1][1].axhline(trim_mean(std, proportiontocut=0.1), ls='--', lw=1, c='k', label=f'$\sigma_{{avg}} = ${trim_mean(std, proportiontocut=0.1):.4f}')

	if title != None :
		fig[0].text(0,0, title, va='baseline', ha='left', size='x-small')

	if last_layer :
		fig[1][0].legend(frameon=True, fontsize='small')
		fig[1][1].legend(frameon=True, fontsize='small')
		plt.subplots_adjust(wspace=0, hspace=0.05)
		fig[1][0].axhline(y=0, alpha=0.3, c='k')
		fig[1][0].set_ylabel('BIAS', loc='top')
		fig[1][1].set_xlabel(xlabel, loc='right')
		fig[1][1].set_ylabel('STD', loc='top')
		if saveas != None :
			sfigs.save_figure(fig[0], outpath, saveas, png=True, pdf=False, pkl=False)
		if show :
			plt.show()


def plot_richness_mass(x_, y_, yerr_=None, masscut=0, Nhex=50,
	xlabel=None, ylabel=None, title=None, label=None,
	fit=None, outpath=None, saveas=None, show=False) :

	if yerr_ is None :
		yerr_ = np.ones_like(y_)

	cut = x_ > masscut
	x = x_[cut]
	y = y_[cut]
	yerr = yerr_[cut]

	fig = plt.subplots(1,1, figsize=figsize())
	fig[1].ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	p = fig[1].hexbin(x, y, gridsize=Nhex, norm=mpl.colors.LogNorm(), edgecolors='none')
	fig[1].grid(which='major', axis='both', ls='-', lw=0.5, color='grey')

	if fit != None :
		x4plot = np.linspace(13,16,100)
		
		params, pcov = curve_fit(globals()[fit], x, y, sigma=yerr)
		perrs = np.sqrt(np.diag(pcov))
		
		if fit == 'linear' :
			slope = params[0]
			inter = params[1]
			shift = 0
			slope_err = perrs[0]
			fig[1].plot(x4plot, globals()[fit](x4plot, *params), lw=1, c='r',
				label='$\lambda \propto M_{{200c}}^{{{{{:.2f}}} \pm {{{:.2f}}}}}$'.format(slope, slope_err))
		elif fit == 'lawnchair' :
			slope = params[0]/params[1]
			inter = params[3]
			shift = params[2]
			slope_err = np.sqrt((perrs[0]/params[1])**2 + (params[0]*perrs[1]/params[1]**2)**2)
			fig[1].plot(x4plot, globals()[fit](x4plot, *params), lw=1, c='r');
			fig[1].plot(x4plot, slope*(x4plot - shift) + inter, lw=1, c='r', ls='--',
				label='$\lambda \propto M_{{200c}}^{{{{{:.2f}}} \pm {{{:.2f}}}}}$'.format(slope, slope_err))

	
	if title != None :
		fig[0].text(0,0, title, va='baseline', ha='left', size='x-small')
	if label != None :
	        plt.legend(loc='lower right', frameon=True)
	        
	fig[1].set_xlabel(xlabel, loc='right')
	fig[1].set_ylabel(ylabel, loc='top')
	#fig[1].set_xlim(13, 15.5)
	#fig[1].set_ylim(0.25, 2.5)
	if saveas != None :
		sfigs.save_figure(fig[0], outpath, saveas, png=True, pdf=False, pkl=False)
	if show :
		plt.show()
	return


def plot_richness_mass_zbinned(x_, y_, z_, yerr_=None, masscut=0, zbins=3, Nhex=50,
	xlabel=None, ylabel=None, title=None, label=None,
	fit=None, outpath=None, saveas=None, show=False) :

	if yerr_ is None :
		yerr_ = np.ones_like(y_)

	cut = x_ > masscut
	x = x_[cut]
	y = y_[cut]
	z = z_[cut]
	yerr = yerr_[cut]

	if type(zbins) is int :
		zbins = np.linspace(0, max(z), zbins)
	
	when_zbin = np.digitize(z, zbins) - 1
	Nzbins = max(when_zbin)

	fig = plt.subplots(1, Nzbins, figsize=figsize(Nzbins,1), sharex=True, sharey=True)
	[ax.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2]) for ax in fig[1]]
	[ax.grid(which='major', axis='both', ls='-', lw=0.5, color='grey') for ax in fig[1]]

	if fig != None :
		x4plot = np.linspace(13,16,100)

		if fit == 'linear' :
			params = [curve_fit(linear, x[when_zbin==i], y[when_zbin==i], sigma=yerr[when_zbin==i]) for i in range(Nzbins)]
			perrs = [np.sqrt(np.diag(ps[1])) for ps in params]
			for i in range(Nzbins) :
				slope = params[i][0][0]
				inter = params[i][0][1]
				shift = 0
				slope_err = perrs[i][0]
				fig[1][i].plot(x4plot, linear(x4plot, *params[i][0]), lw=1, c='r',
					label='$\lambda \propto M_{{200c}}^{{{{{:.2f}}} \pm {{{:.2f}}}}}$'.format(slope, slope_err))
				fig[1][i].hexbin(x[when_zbin==i], y[when_zbin==i], gridsize=Nhex, norm=mpl.colors.LogNorm(), edgecolors='none')
				fig[1][i].annotate(f'${{{zbins[i]:.1f}}} < z_{{halo}} \leq {{{zbins[i+1]:.1f}}}$',
					xy=(fig[1][i].axis()[0], fig[1][i].axis()[3]), xycoords='data', va='top', ha='left',
					bbox=dict(boxstyle='round', fc='w', lw=1))
		elif fit == 'lawnchair' :
			params = [curve_fit(lawnchair, x[when_zbin==i], y[when_zbin==i], sigma=yerr[when_zbin==i], bounds=((0.1,0.1,12,0),(2,2,13.5,1.5)))
				for i in range(Nzbins)]
			perrs = [np.sqrt(np.diag(ps[1])) for ps in params]
			for i in range(Nzbins) :
				slope = params[i][0][0]/params[i][0][1]
				inter = params[i][0][3]
				shift = params[i][0][2]
				slope_err = np.sqrt((perrs[i][0]/params[i][0][1])**2 + (params[i][0][0]*perrs[i][1]/params[i][0][1]**2)**2)
				fig[1][i].plot(x4plot, lawnchair(x4plot, *params[i][0]), lw=1, c='r')
				fig[1][i].plot(x4plot, slope*(x4plot - shift) + inter, lw=1, c='r', ls='--',
					label='$\lambda \propto M_{{200c}}^{{{{{:.2f}}} \pm {{{:.2f}}}}}$'.format(slope, slope_err))
				fig[1][i].hexbin(x[when_zbin==i], y[when_zbin==i], gridsize=Nhex, norm=mpl.colors.LogNorm(), edgecolors='none')
				fig[1][i].annotate(f'${{{zbins[i]:.1f}}} < z_{{halo}} \leq {{{zbins[i+1]:.1f}}}$',
					xy=(fig[1][i].axis()[0], fig[1][i].axis()[3]), xycoords='data', va='top', ha='left',
					bbox=dict(boxstyle='round', fc='w', lw=1))
				fig[1][i].axhline(params[i][0][3], lw=1, ls=':', c='r')
				fig[1][i].axvline(params[i][0][2], lw=1, ls=':', c='r')

	if title != None :
		fig[0].text(0,0, title, va='baseline', ha='left', size='x-small')
	if label != None :
	        [ax.legend(loc='lower right', frameon=True) for ax in fig[1]]
	        
	[ax.set_xlabel(xlabel, loc='right') for ax in fig[1]]
	fig[1][0].set_ylabel(ylabel, loc='top')
	[ax.set_xlim(13, 15.5) for ax in fig[1]]
	[ax.set_ylim(0, 2.5) for ax in fig[1]]
	if saveas != None :
		sfigs.save_figure(fig[0], outpath, saveas, png=True, pdf=False, pkl=False)
	if show :
		plt.show()
	return


def plot_richnessLogProb_binned(rvals, redshift, mass, cfg, xlabel='$\log \lambda$', ylabel='$P\left(\log r | M_{200c}, z\\right)$', title=None,
	outpath=None, saveas=None, show=False) :

	Nzbins = cfg['mass_richness']['binned']['Nzbins']
	NMbins = cfg['mass_richness']['binned']['NMbins']
	Nrbins = cfg['mass_richness']['binned']['Nrbins']
	zbins = np.linspace(cfg['mass_richness']['binned']['zrange'][0], cfg['mass_richness']['binned']['zrange'][1], Nzbins+1)
	Mbins = np.linspace(cfg['mass_richness']['binned']['Mrange'][0], cfg['mass_richness']['binned']['Mrange'][1], NMbins+1)
	rbins = np.linspace(cfg['mass_richness']['binned']['rrange'][0], cfg['mass_richness']['binned']['rrange'][1], Nrbins+1)
	
	where_z_isin = np.digitize(redshift, bins=zbins, right=True) - 1
	where_M_isin = np.digitize(mass, bins=Mbins, right=True) - 1
	
	fig, axs = plt.subplots(NMbins, Nzbins, figsize=figsize(Nzbins-1, NMbins-1))
	fig.subplots_adjust(wspace=0.05, hspace=0.05, left=0.3, bottom=0.3)

	bin_fits = {'Mbin':[], 'zbin':[], 'mu':[], 'mu_err':[], 'sigma':[], 'sigma_err':[]}

	for Mbin in range(NMbins) :
		for zbin in range(Nzbins) :
			ax = axs[Mbin][zbin]
			vals_inMzbin = rvals[(where_z_isin==zbin) & (where_M_isin==Mbin)]
			
			h = ax.hist(vals_inMzbin, bins=rbins, histtype='step', density=True, color='k')
			
			if sum(vals_inMzbin) > cfg['mass_richness']['binned']['min_Ncl_inbin'] :
				try :
					rs = (rbins[1:] + rbins[:-1]) / 2
					
					popt, pcov = curve_fit(gauss, xdata=rs, ydata=h[0])
					a, mu, sigma = popt
					da, dmu, dsigma = np.sqrt(np.diag(pcov))
					
					bin_fits['Mbin'].append([(Mbins[1:] + Mbins[:-1])[Mbin]/2, np.diff(Mbins)[Mbin]/2])
					bin_fits['zbin'].append([(zbins[1:] + zbins[:-1])[zbin]/2, np.diff(zbins)[zbin]/2])
					bin_fits['mu'].append(mu)
					bin_fits['mu_err'].append(dmu)
					bin_fits['sigma'].append(sigma)
					bin_fits['sigma_err'].append(dsigma)

					ax.plot(rs, gauss(rs, *popt), c='C0', lw=2, alpha=0.7)
					ax.axvline(mu, lw=1, c='C1', label=f"$\mu = {{{mu:.2f}}}$")
					ax.fill_between(rs[(rs >= (mu - sigma)) & (rs <= (mu + sigma))],
						gauss(rs[(rs >= (mu - sigma)) & (rs <= (mu + sigma))], *popt),
						color='C2', alpha=0.8, label=f"$\sigma = {{{popt[2]:.2f}}}$")
					leg = ax.legend(labelcolor='mec', frameon=True, fontsize='x-small', loc='best', handlelength=0, handletextpad=0)
					for item in leg.legendHandles :
						item.set_visible(False)
				except :
					pass
			
			ax.set_ylim(bottom=0)
			ax.set_yticklabels([])
			ax.set_yticks([])
			
			if (Mbin == 0) and (zbin == 0) :
				l, b, w, h = ax.get_position().bounds
				plt.figtext(l, b+h, ylabel, va='top', ha='right', rotation='vertical')
			elif (Mbin == NMbins-1) and (zbin == Nzbins-1) :
				l, b, w, h = ax.get_position().bounds
				plt.figtext(l+w, b, xlabel, va='top', ha='right');
			
			if zbin == Nzbins-1 :
				ax.set_ylabel(f"$({{{Mbins[Mbin]:.1f}}},{{{Mbins[Mbin+1]:.1f}}}]$", loc='center', fontsize='small', rotation=270, labelpad=15)
				ax.yaxis.set_label_position('right')
			if Mbin == 0 :
				ax.set_xlabel(f"$({{{zbins[zbin]:.1f}}}, {{{zbins[zbin+1]:.1f}}}]$", loc='right')
				ax.xaxis.set_label_position('top')
	
	if title != None :
		fig.text(0,0, title, va='baseline', ha='left', size='x-small')
	
	if (outpath != None) & (saveas != None) :
		sfigs.save_figure(fig, outpath, saveas=saveas)
	if show :
		plt.show()

	return bin_fits

	
	

def plot_richnessLogProb_binned_popt(bin_fits, bintype='zbin', ylabel='Mo[$P\left(\ln r | M_{200c},z\\right)$]', title=None,
	outpath=None, saveas=None, show=False) :

	if bintype == 'zbin' :
		xlabel = '$\log M_{200c}$'
		leglabel = '$z$'
		xtype = 'Mbin'
	elif bintype == 'Mbin' :
		xlabel = '$\log \lambda$'
		leglabel = '$\log M_{200c}$'
		xtype = 'zbin'

	fig, axs = plt.subplots(1,1, figsize=figsize())

	## collect in either zbin or Mbin
	for ibin, index in zip(*np.unique(np.vstack(bin_fits[bintype])[:,0], return_index=True)) :
		in_bin = (np.vstack(bin_fits[bintype])[:,0] == ibin)

		x = np.vstack(bin_fits[xtype])[:,0][in_bin]
		xerr = np.vstack(bin_fits[xtype])[:,1][in_bin]

		y = np.array(bin_fits['mu'])[in_bin]
		yerr = np.array(bin_fits['sigma'])[in_bin]

		label = f"$({bin_fits[bintype][index][0]-bin_fits[bintype][index][1]:.2f}, {bin_fits[bintype][index][0]+bin_fits[bintype][index][1]:.2f}]$"

		axs.errorbar(x, y, xerr=xerr, yerr=yerr, marker='o', ms=3, label=label, capsize=0.5, capthick=1, elinewidth=1)

	axs.grid(which='major', axis='both', ls='-', lw=0.5, color='grey')

	leg = axs.legend(frameon=True, fontsize='x-small', loc='best', handlelength=0, handletextpad=1, markerscale=0.5, title=leglabel, alignment='right')
	for item, text in zip(leg.legendHandles, leg.get_texts()) :
		item.set_visible(False)
		col = item.get_color()
		if isinstance(col, np.ndarray) :
			col = col[0]
		text.set_color(col)

	axs.set_xlabel(xlabel, loc='right')
	axs.set_ylabel(ylabel, loc='top')

	if title != None :
		fig.text(0,0, title, va='baseline', ha='left', size='x-small')
	
	if (outpath != None) & (saveas != None) :
		sfigs.save_figure(fig, outpath, saveas=saveas)
	if show :
		plt.show()



def plot_richness_richness(x_, y_, yerr_=None, richcut=0, Nhex=50, 
	xlabel=None, ylabel=None, title=None, label=None,
	fit=None, outpath=None, saveas=None, show=False,) :

	if yerr_ is None :
		yerr_ = np.ones_like(y_)
	
	cut = x_ > richcut
	x = x_[cut]
	y = y_[cut]
	yerr = yerr_[cut]
	
	fig = plt.subplots(1,1, figsize=figsize())
	fig[1].ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	p = fig[1].hexbin(x, y, gridsize=Nhex, norm=mpl.colors.LogNorm(), edgecolors='none');
	fig[1].grid(which='major', axis='both', ls='-', lw=0.5, color='grey')

	#if fit != None :
#		params, pcov = curve_fit(linear, x, y, sigma=yerr)
#		perrs = np.sqrt(np.diag(pcov))
#		x4plot = np.linspace(0,3,3)
#		fig[1].plot(x4plot, linear(x4plot, *params), lw=1, c='r',
#			label='$\lambda_|=\lambda_{{-}}^{{({:.2f}\pm{:.2f})}}10^{{({:.2f}\pm{:.2f})}} $'.format(params[1],perrs[1],params[0], perrs[0]));
	if fit != None :
		x4plot = np.linspace(0,3,3)
		
		params, pcov = curve_fit(linear, x, y)#, sigma=yerr, bounds=((0,-0.5),(2,0.5)))
		perrs = np.sqrt(np.diag(pcov))

		if title != None :
			cats = title.split(' - ')
			cat1 = cats[0].split('.')[0]
			cat2 = cats[1].split('.')[0]
		else :
			cat1 = '_'
			cat2 = '|'
		
		slope = params[0]
		inter = params[1]
		slope_err = perrs[0]
		inter_err = perrs[1]
		label = f"$\lambda_{{{cat2}}} = 10^{{{inter:.2f}}}\lambda_{{{cat1}}}^{{{slope:.2f}}}$"
		fig[1].plot(x4plot, linear(x4plot, *params), lw=1, c='r', label=label)

	
	if title != None :
		fig[0].text(0,0, title, va='baseline', ha='left', size='x-small')
	if label != None :
	        plt.legend(loc='lower right', frameon=True)
	        
	fig[1].set_xlabel(xlabel, loc='right')
	fig[1].set_ylabel(ylabel, loc='top')
	fig[1].set_xlim(0, 3)
	fig[1].set_ylim(0,3)


	if (outpath != None) & (saveas != None) :
		sfigs.save_figure(fig[0], outpath, saveas=saveas)
	if show :
		plt.show()



def plot_richness_richness_zbinned(x_, y_, z_, yerr_=None, richcut=0, zbins=3, Nhex=50,
	xlabel=None, ylabel=None, title=None, label=None,
	fit=None, outpath=None, saveas=None, show=False) :

	if yerr_ is None :
		yerr_ = np.ones_like(y_)

	cut = x_ > richcut
	x = x_[cut]
	y = y_[cut]
	z = z_[cut]
	yerr = yerr_[cut]

	if type(zbins) is int :
		zbins = np.linspace(0, max(z), zbins)
	
	when_zbin = np.digitize(z, zbins) - 1
	Nzbins = max(when_zbin)

	fig = plt.subplots(1, Nzbins, figsize=figsize(Nzbins,1), sharex=True, sharey=True)
	[ax.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2]) for ax in fig[1]]
	[ax.grid(which='major', axis='both', ls='-', lw=0.5, color='grey') for ax in fig[1]]

	if fig != None :
		x4plot = np.linspace(0,3,3)

		if fit != None :
			if title != None :
				cats = title.split(' - ')
				cat1 = cats[0].split('.')[0]
				cat2 = cats[1].split('.')[0]
			else :
				cat1 = '_'
				cat2 = '|'
		
			params = [curve_fit(linear, x[when_zbin==i], y[when_zbin==i]) for i in range(Nzbins)]
			perrs = [np.sqrt(np.diag(ps[1])) for ps in params]
			for i in range(Nzbins) :
				slope = params[i][0][0]
				inter = params[i][0][1]
				shift = 0
				slope_err = perrs[i][0]
				label = f"$\lambda_{{{cat2}}} = 10^{{{inter:.2f}}}\lambda_{{{cat1}}}^{{{slope:.2f}}}$"
				fig[1][i].plot(x4plot, linear(x4plot, *params[i][0]), lw=1, c='r',label=label)
				fig[1][i].hexbin(x[when_zbin==i], y[when_zbin==i], gridsize=Nhex, norm=mpl.colors.LogNorm(), edgecolors='none')
				fig[1][i].annotate(f'${{{zbins[i]:.1f}}} < z_{{halo}} \leq {{{zbins[i+1]:.1f}}}$',
					xy=(fig[1][i].axis()[0], fig[1][i].axis()[3]), xycoords='data', va='top', ha='left',
					bbox=dict(boxstyle='round', fc='w', lw=1))

	#if title != None :
	#	fig[0].text(0,0, title, va='baseline', ha='left', size='x-small')
	if label != None :
	        [ax.legend(loc='lower right', frameon=True) for ax in fig[1]]
	        
	[ax.set_xlabel(xlabel, loc='right') for ax in fig[1]]
	fig[1][0].set_ylabel(ylabel, loc='top')
	[ax.set_xlim(0, 3) for ax in fig[1]]
	[ax.set_ylim(0, 3) for ax in fig[1]]
	if saveas != None :
		sfigs.save_figure(fig[0], outpath, saveas, png=True, pdf=False, pkl=False)
	if show :
		plt.show()
