
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


def plot_richness_mass(x_, y_, yerr_=None, masscut=0, Nhex=50,
	xlabel=None, ylabel=None, title=None, fit=None,
	outpath='./', saveas='richness_mass', show=False) :
	'''
	Produces a hexbin plot of the richness-mass distribution.
	
	Attributes
	----------
	(x_, y_) : array_like
		The mass (x_) and richness (y_).
	yerr_ : array_like
		The richness error.
	masscut : float
		A mass cut can be applied to ignore halos with mass below a certain threshold.
	Nhex : int
		The granularity of the hexbins.
	(xlabel, ylabel) : str
		The x- and y-labels for the plot.
	title : str
		Additional descriptor to be placed in lower left corner of the plot.
	fit : str
		The fitting function to be used (linear or lawnchair). If None then no fit will be made.
	outpath : str
		The output path for the plots to be stored. The default is to save in the present work dir.
	saveas : str
		The file name to save the plot under. Default is richness_mass.
	show : bool
		If the plot should be displayed upon running.	
	'''

	if yerr_ is None :
		yerr_ = np.ones_like(y_)

	cut = x_ > masscut
	x = x_[cut]
	y = y_[cut]
	yerr = yerr_[cut]

	fig, axs = plt.subplots(1,1, figsize=figsize())
	axs.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	p = axs.hexbin(x, y, gridsize=Nhex, norm=mpl.colors.LogNorm(), edgecolors='none')
	axs.grid(which='major', axis='both', ls='-', lw=0.5, color='grey')

	if fit != None :
		x4plot = np.linspace(13,16,100)
		
		params, pcov = curve_fit(globals()[fit], x, y, sigma=yerr)
		perrs = np.sqrt(np.diag(pcov))
		
		if fit == 'linear' :
			slope = params[0]
			inter = params[1]
			shift = 0
			slope_err = perrs[0]

			label='$\lambda \propto M_{{200c}}^{{{{{:.2f}}} \pm {{{:.2f}}}}}$'.format(slope, slope_err)
			axs.plot(x4plot, globals()[fit](x4plot, *params), lw=1, c='r', label=label)
		elif fit == 'lawnchair' :
			slope = params[0]/params[1]
			inter = params[3]
			shift = params[2]
			slope_err = np.sqrt((perrs[0]/params[1])**2 + (params[0]*perrs[1]/params[1]**2)**2)

			axs.plot(x4plot, globals()[fit](x4plot, *params), lw=1, c='r');
			label='$\lambda \propto M_{{200c}}^{{{{{:.2f}}} \pm {{{:.2f}}}}}$'.format(slope, slope_err)
			axs.plot(x4plot, slope*(x4plot - shift) + inter, lw=1, c='r', ls='--', label=label)

	
	if title != None :
		fig.text(0,0, title, va='baseline', ha='left', size='x-small')
	if fit != None :
	        plt.legend(loc='lower right', frameon=True)
	        
	axs.set_xlabel(xlabel, loc='right')
	axs.set_ylabel(ylabel, loc='top')
	#axs.set_xlim(13, 15.5)
	#axs.set_ylim(0.25, 2.5)

	if saveas != None :
		sfigs.save_figure(fig, outpath, saveas)
	if show :
		plt.show()
	return


def plot_richness_mass_zbinned(x_, y_, z_, yerr_=None, masscut=0, zbins=3, Nhex=50,
	xlabel=None, ylabel=None, title=None, fit=None, outpath=None, saveas=None, show=False) :
	'''
	Produces a series of hexbin plots of the richness-mass distribution binned by redshift.
	
	Attributes
	----------
	(x_, y_) : array_like
		The mass (x_) and richness (y_).
	yerr_ : array_like
		The richness error.
	z_ : array_like
		Redshifts
	masscut : float
		A mass cut can be applied to ignore halos with mass below a certain threshold.
	zbins : int
		The number of redshift bins to use. This number of subplots will be made.
	Nhex : int
		The granularity of the hexbins.
	(xlabel, ylabel) : str
		The x- and y-labels for the plot.
	title : str
		Additional descriptor to be placed in lower left corner of the plot.
	fit : str
		The fitting function to be used (linear or lawnchair). If None then no fit will be made.
	outpath : str
		The output path for the plots to be stored. The default is to save in the present work dir.
	saveas : str
		The file name to save the plot under. Default is richness_mass.
	show : bool
		If the plot should be displayed upon running.
	'''

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

	fig, axs = plt.subplots(1, Nzbins, figsize=figsize(Nzbins,1), sharex=True, sharey=True)
	[ax.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2]) for ax in axs]
	[ax.grid(which='major', axis='both', ls='-', lw=0.5, color='grey') for ax in axs]

	if fig != None :
		x4plot = np.linspace(13,16,100)

		if fit == 'linear' :
			params = [curve_fit(linear, x[when_zbin==i], y[when_zbin==i], sigma=yerr[when_zbin==i])
					for i in range(Nzbins)]
			perrs = [np.sqrt(np.diag(ps[1])) for ps in params]
			for i in range(Nzbins) :
				slope = params[i][0][0]
				inter = params[i][0][1]
				shift = 0
				slope_err = perrs[i][0]

				label = '$\lambda \propto M_{{200c}}^{{{{{:.2f}}} \pm {{{:.2f}}}}}$'.format(slope, slope_err)
				axs[i].plot(x4plot, linear(x4plot, *params[i][0]), lw=1, c='r', label=label)
				axs[i].hexbin(x[when_zbin==i], y[when_zbin==i], gridsize=Nhex, norm=mpl.colors.LogNorm(), edgecolors='none')
				axs[i].annotate(f'${{{zbins[i]:.1f}}} < z_{{halo}} \leq {{{zbins[i+1]:.1f}}}$',
					xy=(axs[i].axis()[0], axs[i].axis()[3]), xycoords='data', va='top', ha='left',
					bbox=dict(boxstyle='round', fc='w', lw=1))
		elif fit == 'lawnchair' :
			bounds = ((0.1,0.1,12,0),(2,2,13.5,1.5))
			params = [curve_fit(lawnchair, x[when_zbin==i], y[when_zbin==i], sigma=yerr[when_zbin==i], bounds=bounds)
				for i in range(Nzbins)]
			perrs = [np.sqrt(np.diag(ps[1])) for ps in params]
			for i in range(Nzbins) :
				slope = params[i][0][0]/params[i][0][1]
				inter = params[i][0][3]
				shift = params[i][0][2]
				slope_err = np.sqrt((perrs[i][0]/params[i][0][1])**2
						+ (params[i][0][0]*perrs[i][1]/params[i][0][1]**2)**2)

				axs[i].plot(x4plot, lawnchair(x4plot, *params[i][0]), lw=1, c='r')
				label='$\lambda \propto M_{{200c}}^{{{{{:.2f}}} \pm {{{:.2f}}}}}$'.format(slope, slope_err)
				axs[i].plot(x4plot, slope*(x4plot - shift) + inter, lw=1, c='r', ls='--', label=label)
				axs[i].hexbin(x[when_zbin==i], y[when_zbin==i], gridsize=Nhex, norm=mpl.colors.LogNorm(), edgecolors='none')
				axs[i].annotate(f'${{{zbins[i]:.1f}}} < z_{{halo}} \leq {{{zbins[i+1]:.1f}}}$',
					xy=(axs[i].axis()[0], axs[i].axis()[3]), xycoords='data', va='top', ha='left',
					bbox=dict(boxstyle='round', fc='w', lw=1), zorder=np.inf)
				axs[i].axhline(params[i][0][3], lw=1, ls=':', c='r')
				axs[i].axvline(params[i][0][2], lw=1, ls=':', c='r')

	if title != None :
		fig.text(0,0, title, va='baseline', ha='left', size='x-small')
	if label != None :
	        [ax.legend(loc='lower right', frameon=True) for ax in axs]
	        
	[ax.set_xlabel(xlabel, loc='right') for ax in axs]
	axs[0].set_ylabel(ylabel, loc='top')
	[ax.set_xlim(13, 15.5) for ax in axs]
	[ax.set_ylim(0, 2.5) for ax in axs]

	if saveas != None :
		sfigs.save_figure(fig, outpath, saveas)
	if show :
		plt.show()
	return


def plot_richnessLogProb_binned(rvals, redshift, mass, cfg, xlabel='$\log \lambda$', ylabel='$P\left(\log r | M_{200c}, z\\right)$',
	title=None, outpath=None, saveas=None, show=False) :
	'''
	Produces an array of fitted histograms for the different mass and redshift bins.
	
	Attributes
	----------
	rvals : array_like
		Richness
	redshift : array_like
		Redshift
	mass : array_like
		Halo masses
	cfg : dict
		The configuration dictionary.
	(xlabel, ylabel) : str
		The x- and y-labels for the plot.
	title : str
		Additional descriptor to be placed in lower left corner of the plot.
	outpath : str
		The output path for the plots to be stored. The default is to save in the present work dir.
	saveas : str
		The file name to save the plot under. Default is richness_mass.
	show : bool
		If the plot should be displayed upon running.
	'''


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
				label = f"$({{{Mbins[Mbin]:.1f}}},{{{Mbins[Mbin+1]:.1f}}}]$"
				ax.set_ylabel(label, loc='center', fontsize='small', rotation=270, labelpad=15)
				ax.yaxis.set_label_position('right')
			if Mbin == 0 :
				label = f"$({{{zbins[zbin]:.1f}}}, {{{zbins[zbin+1]:.1f}}}]$"
				ax.set_xlabel(label, loc='right')
				ax.xaxis.set_label_position('top')
	
	if title != None :
		fig.text(0,0, title, va='baseline', ha='left', size='x-small')
	
	if (outpath != None) & (saveas != None) :
		sfigs.save_figure(fig, outpath, saveas=saveas)
	if show :
		plt.show()

	return bin_fits

	
	

def plot_richnessLogProb_binned_popt(bin_fits, bintype='zbin', ylabel='MODE[$P\left(\ln r | M_{200c},z\\right)$]', title=None,
	outpath=None, saveas=None, show=False) :
	'''
	Plots the best fit parameters from the binned analysis of the richness-mass plane.
	
	Attributes
	----------
	rvals : array_like
		Richness
	redshift : array_like
		Redshift
	mass : array_like
		Halo masses
	cfg : dict
		The configuration dictionary.
	(xlabel, ylabel) : str
		The x- and y-labels for the plot.
	title : str
		Additional descriptor to be placed in lower left corner of the plot.
	outpath : str
		The output path for the plots to be stored. The default is to save in the present work dir.
	saveas : str
		The file name to save the plot under. Default is richness_mass.
	show : bool
		If the plot should be displayed upon running.	
	'''

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



def plot_richness_richness(r1_, r2_, r1err_=None, r2err_=None, richcut1=0, richcut2=0, Nhex=50, 
	xlabel=None, ylabel=None, title=None, fit=None, outpath=None, saveas=None, show=False,) :
	'''
	Produces a hexbin plot of the richness-richness space.
	
	Attributes
	----------
	(r1_, r2_) : array_like
		Richness for cat1 and cat2
	(r1err_, r2err_) : array_like
		Richness errors for cat1 and cat2
	(richcut1, richcut2) : float
		The richness cuts to have applied to the catalogs.
	Nhex : int
		The granularity of the hexbins.
	(xlabel, ylabel) : str
		The x- and y-labels for the plot.
	title : str
		Additional descriptor to be placed in lower left corner of the plot.
	fit : str
		The fitting function to use -- see fit_func.py. If None then no fit is produced.
	outpath : str
		The output path for the plots to be stored. The default is to save in the present work dir.
	saveas : str
		The file name to save the plot under. Default is richness_mass.
	show : bool
		If the plot should be displayed upon running.	
	'''

	if r1err_ is None :
		r1err_ = np.ones_like(r1_)
	if r2err_ is None :
		r2err_ = np.ones_like(r2_)
	
	cut = (r1_ > richcut1) & (r2_ > richcut2)
	r1 = r1_[cut]
	r2 = r2_[cut]
	r1err = r1err_[cut]
	r2err = r2err_[cut]
	
	fig, axs = plt.subplots(1,1, figsize=figsize())
	axs.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	p = axs.hexbin(r1, r2, gridsize=Nhex, norm=mpl.colors.LogNorm(), edgecolors='none');
	axs.grid(which='major', axis='both', ls='-', lw=0.5, color='grey')

	if fit != None :
		x4plot = np.linspace(0,3,3)
		
		params, pcov = curve_fit(linear, r1, r2)#, sigma=yerr, bounds=((0,-0.5),(2,0.5)))
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
		axs.plot(x4plot, linear(x4plot, *params), lw=1, c='r', label=label)

	
	if title != None :
		fig.text(0,0, title, va='baseline', ha='left', size='x-small')
	if label != None :
	        plt.legend(loc='lower right', frameon=True)
	        
	axs.set_xlabel(xlabel, loc='right')
	axs.set_ylabel(ylabel, loc='top')
	axs.set_xlim(0, 3)
	axs.set_ylim(0,3)


	if (outpath != None) & (saveas != None) :
		sfigs.save_figure(fig, outpath, saveas=saveas)
	if show :
		plt.show()



def plot_richness_richness_zbinned(r1_, r2_, z_, r1err_=None, r2err_=None, richcut1=0, richcut2=0, zbins=3, Nhex=50,
	xlabel=None, ylabel=None, title=None, fit=None, outpath=None, saveas=None, show=False) :
	'''
	Produces a series of hexbin plots of the richness-richness space binned in redshift.
	
	Attributes
	----------
	(r1_, r2_) : array_like
		Richness for cat1 and cat2
	z_ : array_like
		Redshift
	(r1err_, r2err_) : array_like
		Richness errors for cat1 and cat2
	(richcut1, richcut2) : float
		The richness cuts to have applied to the catalogs.
	zbins : int
		The number of redshift bins to use. This number of subplots will be made.
	Nhex : int
		The granularity of the hexbins.
	(xlabel, ylabel) : str
		The x- and y-labels for the plot.
	title : str
		Additional descriptor to be placed in lower left corner of the plot.
	fit : str
		The fitting function to use -- see fit_func.py. If None then no fit is produced.
	outpath : str
		The output path for the plots to be stored. The default is to save in the present work dir.
	saveas : str
		The file name to save the plot under. Default is richness_mass.
	show : bool
		If the plot should be displayed upon running.	
	'''

	if r1err_ is None :
		r1err_ = np.ones_like(r1_)
	if r2err_ is None :
		r2err_ = np.ones_like(r2_)

	cut = (r1_ > richcut1) & (r2_ > richcut2)
	r1 = r1_[cut]
	r2 = r2_[cut]
	z = z_[cut]
	r1err = r1err_[cut]
	r2err = r2err_[cut]

	if type(zbins) is int :
		zbins = np.linspace(0, max(z), zbins)
	
	when_zbin = np.digitize(z, zbins) - 1
	Nzbins = max(when_zbin)

	fig, axs = plt.subplots(1, Nzbins, figsize=figsize(Nzbins,1), sharex=True, sharey=True)
	[ax.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2]) for ax in axs]
	[ax.grid(which='major', axis='both', ls='-', lw=0.5, color='grey') for ax in axs]

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
		
			params = [curve_fit(linear, r1[when_zbin==i], r2[when_zbin==i]) for i in range(Nzbins)]
			perrs = [np.sqrt(np.diag(ps[1])) for ps in params]
			for i in range(Nzbins) :
				slope = params[i][0][0]
				inter = params[i][0][1]
				shift = 0
				slope_err = perrs[i][0]
				label = f"$\lambda_{{{cat2}}} = 10^{{{inter:.2f}}}\lambda_{{{cat1}}}^{{{slope:.2f}}}$"
				axs[i].plot(x4plot, linear(x4plot, *params[i][0]), lw=1, c='r',label=label)
				axs[i].hexbin(r1[when_zbin==i], r2[when_zbin==i], gridsize=Nhex, norm=mpl.colors.LogNorm(), edgecolors='none')
				axs[i].annotate(f'${{{zbins[i]:.1f}}} < z_{{halo}} \leq {{{zbins[i+1]:.1f}}}$',
					xy=(axs[i].axis()[0], axs[i].axis()[3]), xycoords='data', va='top', ha='left',
					bbox=dict(boxstyle='round', fc='w', lw=1), zorder=np.inf)

	#if title != None :
	#	fig[0].text(0,0, title, va='baseline', ha='left', size='x-small')
	if label != None :
	        [ax.legend(loc='lower right', frameon=True) for ax in axs]
	        
	[ax.set_xlabel(xlabel, loc='right') for ax in axs]
	axs[0].set_ylabel(ylabel, loc='top')
	[ax.set_xlim(0, 3) for ax in axs]
	[ax.set_ylim(0, 3) for ax in axs]
	if saveas != None :
		sfigs.save_figure(fig, outpath, saveas)
	if show :
		plt.show()
