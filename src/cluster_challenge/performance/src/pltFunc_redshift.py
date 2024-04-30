
import numpy as np

from scipy.stats import binned_statistic, trim_mean
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

from plot_style import set_style, figsize
import saving_figures as sfigs
from fit_func import linear, lawnchair, gauss


## Initialize style for plots.
set_style()




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


