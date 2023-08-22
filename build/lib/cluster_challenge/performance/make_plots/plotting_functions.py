
import numpy as np
import glob
import os
import struct

from scipy.stats import binned_statistic, trim_mean
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.optimize import minimize, curve_fit
from scipy.special import erf


def save_figure(fig, outpath, parent_dir, saveas) :
	for path in outpath :
		if not os.path.exists(path+parent_dir) :
			print(f'Making directory:\t{path}{parent_dir}')
			os.makedirs(path+parent_dir)
	fig.savefig(outpath[0] + parent_dir + saveas + '.png', bbox_inches='tight')
	fig.savefig(outpath[1] + parent_dir + saveas + '.pdf', bbox_inches='tight')
	plt.close(fig)


def update_index_file_for_html_display(outpath, category, description) :
	files = list(filter(os.path.isfile, glob.glob(f'{outpath}png_plots/{category}/*')))
	files.sort(key=lambda x: os.path.getmtime(x))

	width = []
	height = []
	for f in files :
		with open(f, 'rb') as fhandle :
			head = fhandle.read(24)
			w, h = struct.unpack('>ii', head[16:24])
			width.append(w)
			height.append(h)

	if not os.path.exists(f'{outpath}png_plots/{category}/display/') :
		print(f'Making directory:\t{outpath}png_plots/{category}/display/')
		os.makedirs(f'{outpath}png_plots/{category}/display/')

	files = [f"https://me.lsst.eu/rsolomon{file.removeprefix('/sps/lsst/users/rsolomon/web')}" for file in files]

	text  = f'<html>\n'
	text += f'<h2>{description}</h2>\n\n'

	text += f'<ul>\n'
	text += '\n'.join([f'<IMG align=center width={w} height={h} SRC={SRC}>' for w,h,SRC in zip(width,height,files)])
	text += f'\n</ul>\n\n'

	text += f'</html>'

	with open(f'{outpath}png_plots/{category}/display/index.html', 'w') as f :
		f.write(text)
	


## DEFINE FITTING FUNCTION
def fit_linear(B, x, y, error) :
	ypred = B[0] + B[1]*x
	s2 = sum([d**2 for d in error])
	negLL = -np.sum(-(y - ypred)**2 / 2 / s2 - np.log(np.sqrt(2 * np.pi * s2)))
	return negLL

def fit_quadratic(B, x, y, error) :
	ypred = B[0] + B[1]*x + B[2]*x**2
	s2 = sum([d**2 for d in error])
	negLL = -np.sum(-(y - ypred)**2 / 2 / s2 - np.log(np.sqrt(2 * np.pi * s2)))
	return negLL



## DEFINE HEX SIZES FOR HEXBIN
hex_h = 1
hex_w = np.sqrt(3)*hex_h

#results = minimize(MLE, [0,-0.5, 0.5], args=(np.array(m), np.array(logsf), np.array(logsf_errs)), method='Nelder-Mead', bounds=((-1,1),(-5,5),(0,1)))
#print(results.x)



def plot_on_sky(ra, dec, size=0.3, alpha=0.05, color='C0', title=None, label=None, outpath=None, saveas=None, show=False,
	only_layer=True, axs=None, last_layer=True) :
	
	if only_layer :         ## IN CASE OF MULTIPLE DATASETS, PLACE fig, axs = ... BEFORE FUNCTION CALL
	        fig, axs = plt.subplots(1,1, figsize=(3.5,3.5))
	
	#axs.scatter(ra, dec, s=size, alpha=alpha, c=color, label=label)
	Nhex = 200
	axs.hexbin(ra, dec, gridsize=(int(Nhex*hex_w),int(Nhex*hex_h)), norm=mpl.colors.LogNorm());
	
	if title != None :
		fig.text(0,0, title, va='baseline', ha='left', size='small')
	if last_layer and label != None :
		axs.annotate(label, xy=(axs.axis()[1], axs.axis()[3]), xycoords='data', va='top', ha='right', bbox=dict(boxstyle='round', fc='w', lw=1))
	if last_layer :
		axs.set_xlabel('RA', loc='right')
		axs.set_ylabel('DEC', loc='top')
		if saveas != None :
			save_figure(fig, outpath, 'onsky/', saveas)
			#for path in outpath :
			#	if not os.path.exists(path+'onsky/') :
			#		os.makedirs(path+'onsky/')
			#
			#plt.savefig(outpath[0] + 'onsky/' + saveas + '.png')
			#plt.savefig(outpath[1] + 'onsky/' + saveas + '.pdf')
		if show :
			plt.show()



def plot_redshift_hist(zs, bin_width=0.05, fill=False, color='C0', title=None, label=None, annotate=None, outpath=None, saveas=None, show=False,
	only_layer=True, fig=None, last_layer=True) :
	
	#if only_layer :        ## IN CASE OF MULTIPLE DATASETS ON SAME PLOT
	if fig == None :        ## IN CASE OF MULTIPLE DATASETS ON SAME PLOT
		fig = plt.subplots(1,1, figsize=(3.5,3.5))
	fig[1].ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	if fill == True :
	        alpha = 0.3
	else :
	        alpha = 1

	bins = int(np.ptp(zs)/bin_width)
	h = fig[1].hist(zs, bins=bins, histtype='step', fill=fill, color=color, alpha=alpha, label=label)
	
	#if title != None :
	#	fig.text(0,0, title, va='baseline', ha='left')
	if last_layer and label != None :
		fig[1].set_ylim(top=1.3*np.max(fig[1].get_ylim()))
		leg = fig[1].legend(labelcolor='mec', frameon=False, fontsize='small', loc='upper left', handlelength=0, handletextpad=0)
		for item in leg.legendHandles :
			item.set_visible(False)
	if last_layer and annotate != None :
		fig[1].annotate(annotate, xy=(fig[1].axis()[1], fig[1].axis()[3]), xycoords='data', va='top', ha='right',
			bbox=dict(boxstyle='round', fc='w', lw=1), zorder=np.inf)
	if last_layer :
		fig[1].set_xlabel('Redshift', loc='right')
		fig[1].set_ylabel('Count', loc='top')
		if title != None :
			fig[0].text(0,0, title, va='baseline', ha='left', size='small')
		if saveas != None :
			save_figure(fig[0], outpath, 'redshifts/', saveas)
		if show :
			plt.show()



def plot_redshift_zVSz(z1, z2, xlabel=None, ylabel=None, diagonal=False, title=None, label=None, outpath=None, saveas=None, show=False,) :
	fig, axs = plt.subplots(1,1, figsize=(3.5,3.5))
	axs.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	if diagonal :
		axs.plot([min(z1),max(z1)],[min(z1),max(z1)], lw=0.5, c='r', ls='--')

	Nhex = 200
	axs.hexbin(z1, z2, gridsize=(int(Nhex*hex_w),int(Nhex*hex_h)), norm=mpl.colors.LogNorm());
	
	if title != None :
		fig.text(0,0, title, va='baseline', ha='left', size='small')
	if label != None :
	        plt.legend()
	        
	axs.set_xlabel(xlabel, loc='right')
	axs.set_ylabel(ylabel, loc='top')
	axs.set_xlim(0,1.5)
	axs.set_ylim(0,1.5)
	if saveas != None :
		save_figure(fig, outpath, 'redshifts/', saveas)
	if show :
		plt.show()


def plot_redshift_std_and_mean(z1, z2, Nbins=50, label=None, title=None, outpath=None, saveas=None, show=False, only_layer=True, fig=None, last_layer=True) :
	zbins = np.linspace(min(z1), max(z1), Nbins)
	vals = (z1-z2)/(1+z1)
	mean = binned_statistic(z1[vals<0.15], vals[vals<0.15], bins=Nbins, statistic='mean')[0]
	std  = binned_statistic(z1[vals<0.15], vals[vals<0.15], bins=Nbins, statistic='std')[0]

	if only_layer :
		fig = plt.subplots(2,1, figsize=(7,3.5), gridspec_kw={'hspace':0.4})
		[ax.ticklabel_format(axis='both', style='scientific', scilimits=[-2,2]) for ax in fig[1]]
	
	fig[1][0].plot(zbins+0.5*np.mean(np.diff(zbins)), mean, lw=1)
	fig[1][0].axhline(trim_mean(mean, proportiontocut=0.1), ls='--', lw=1, c='k', label=f'$b_{{avg}} = ${trim_mean(mean, proportiontocut=0.1):.4f}')
	
	fig[1][1].plot(zbins+0.5*np.mean(np.diff(zbins)), std, lw=1)
	fig[1][1].axhline(trim_mean(std, proportiontocut=0.1), ls='--', lw=1, c='k', label=f'$\sigma_{{avg}} = ${trim_mean(std, proportiontocut=0.1):.4f}')

	if title != None :
		fig[0].text(0,0, title, va='baseline', ha='left', size='small')

	if last_layer :
		fig[1][0].legend(frameon=False, fontsize='small')
		fig[1][1].legend(frameon=False, fontsize='small')
		plt.subplots_adjust(wspace=0, hspace=0.05)
		fig[1][0].axhline(y=0, alpha=0.3, c='k')
		fig[1][0].set_ylabel('BIAS', loc='top')
		fig[1][1].set_xlabel('truez', loc='right')
		fig[1][1].set_ylabel('STD', loc='top')
		if saveas != None :
			save_figure(fig[0], outpath, 'redshifts/', saveas)
			#for path in outpath :
			#	if not os.path.exists(path+'redshifts/') :
			#		os.makedirs(path+'redshifts/')
			#plt.savefig(outpath[0] + 'redshifts/' + saveas + '.png')
			#plt.savefig(outpath[1] + 'redshifts/' + saveas + '.pdf')
		if show :
			plt.show()


def plot_richness_mass(r1, r2, xlabel=None, ylabel=None, title=None, label=None, fit=False, outpath=None, saveas=None, show=False,) :
	fig = plt.subplots(1,1, figsize=(3.5,3.5))
	fig[1].ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	Nhex = 30
	p = fig[1].hexbin(r1, r2, gridsize=(int(Nhex*hex_w),int(Nhex*hex_h)), norm=mpl.colors.LogNorm(), edgecolors='none');

	if fit :
		def richness_mass_relation(x, c0, c1, c2) :
			return 0.5*c0*((x-c1)*erf(x-c1) + 1/np.sqrt(np.pi)*np.exp(-(x-c1)**2) + (x-c1)) + c2

		params, pcov = curve_fit(richness_mass_relation, r1, r2, bounds=((0.5,12,0),(5,15,2)))
		perrs = np.sqrt(np.diag(pcov))
		x4plot = np.linspace(13,16,100)
		fig[1].plot(x4plot, richness_mass_relation(x4plot, *params), lw=1, c='r');
		fig[1].plot(x4plot, params[0]*(x4plot - params[1]) + params[2], lw=1, c='r', ls='--',
			label='$\lambda \propto M_{{200c}}^{{{{{:.2f}}} \pm {{{:.2f}}}}}$'.format(params[0], perrs[0]))

	
	if title != None :
		fig[0].text(0,0, title, va='baseline', ha='left', size='small')
	if label != None :
	        plt.legend(loc='lower right', frameon=False)
	        
	fig[1].set_xlabel(xlabel, loc='right')
	fig[1].set_ylabel(ylabel, loc='top')
	fig[1].set_xlim(13, 15.5)
	fig[1].set_ylim(0.25, 2.5)
	if saveas != None :
		save_figure(fig[0], outpath, 'richness/', saveas)
		#for path in outpath :
		#	if not os.path.exists(path+'richness/') :
		#		os.makedirs(path+'richness/')
		#plt.savefig(outpath[0] + 'richness/' + saveas + '.png')
		#plt.savefig(outpath[1] + 'richness/' + saveas + '.pdf')
	if show :
		plt.show()
	return


def plot_richness_richness(r1, r2, xlabel=None, ylabel=None, title=None, label=None, fit=False, outpath=None, saveas=None, show=False,) :
	fig = plt.subplots(1,1, figsize=(3.5,3.5))
	fig[1].ticklabel_format(axis='both', style='scientific', scilimits=[-2,2])
	
	Nhex = 30
	p = fig[1].hexbin(r1, r2, gridsize=(int(Nhex*hex_w),int(Nhex*hex_h)), norm=mpl.colors.LogNorm(), edgecolors='none');

	if fit :
		def richness_richness_relation(x, c0, c1) :
			return c0 + c1*x

		params, pcov = curve_fit(richness_richness_relation, r1, r2, bounds=((-2,-2),(2,2)))
		perrs = np.sqrt(np.diag(pcov))
		x4plot = np.linspace(0,3,3)
		fig[1].plot(x4plot, richness_richness_relation(x4plot, *params), lw=1, c='r',
			label='$\lambda_|=\lambda_{{-}}^{{({:.2f}\pm{:.2f})}}10^{{({:.2f}\pm{:.2f})}} $'.format(params[1],perrs[1],params[0], perrs[0]));

	
	if title != None :
		fig[0].text(0,0, title, va='baseline', ha='left', size='small')
	if label != None :
	        plt.legend(loc='lower right', frameon=False)
	        
	fig[1].set_xlabel(xlabel, loc='right')
	fig[1].set_ylabel(ylabel, loc='top')
	fig[1].set_xlim(0.5, 3)
	fig[1].set_ylim(0.25,3)
	if saveas != None :
		save_figure(fig[0], outpath, 'richness/', saveas)
		#for path in outpath :
		#	if not os.path.exists(path+'richness/') :
		#		os.makedirs(path+'richness/')
		#plt.savefig(outpath[0] + 'richness/' + saveas + '.png')
		#plt.savefig(outpath[1] + 'richness/' + saveas + '.pdf')
	if show :
		plt.show()
	return
