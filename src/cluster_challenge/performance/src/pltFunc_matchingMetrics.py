
import numpy as np
import math
from src.plot_style import figsize

from clevar.match_metrics.recovery import ClCatalogFuncs as r_cf
from clevar.match_metrics.recovery.array_funcs import get_recovery_rate

import matplotlib.pyplot as plt



def calcTPR(x, prob, cut) :
	y = (x != None)
	num = sum(y & (prob >= cut))
	den = sum(y)
	return num / den

def calcFPR(x, prob, cut) :
	y = (x == None)
	num = sum(y & (prob < cut))
	den = sum(y)
	return 1 - num / den

def make_ROC_curve(cat, steps=np.linspace(1,20,101)) :
	TPR = np.array([calcTPR(cat['mt_cross'], cat['snr_cl'], snr) for snr in steps])
	FPR = np.array([calcFPR(cat['mt_cross'], cat['snr_cl'], snr) for snr in steps])

	fig, axs = plt.subplots(1,1, figsize=figsize())
	axs.scatter(FPR[::10], TPR[::10], s=3)
	axs.plot(FPR, TPR, lw=1)
	axs.plot([0,1],[0,1], lw=1, ls='--', c='k')
	optimalSNR = np.linspace(1,20,101)[np.argmax((1-FPR)*TPR)]
	axs.scatter(FPR[np.argmax((1-FPR)*TPR)], TPR[np.argmax((1-FPR)*TPR)], c='k', label=f"SNR = {optimalSNR:.2f}")
	axs.set_xlim(0,1)
	axs.set_ylim(0,1)
	axs.set_xlabel('FPR', loc='right')
	axs.set_ylabel('TPR', loc='top')
	[axs.text(FPR[i], TPR[i], f"  {np.linspace(1,20,101)[i]:.1f}") for i in range(101)[::10][2:-2]]
	axs.legend(loc='lower right')
	fig.text(0,0, f"AUC = {np.trapz(TPR[::-1], FPR[::-1]):.4f}", va='baseline', ha='left')

	return fig, np.vstack([steps,TPR]).T, np.vstack([steps,FPR]).T


def make_PurityCompleteness_curve(
	hls, cls,
	steps=np.logspace(0,2,201),
	massRange=[1e14,1e17],
	richRange=[0,300],
	speczRange=[0,1.5],
	photozRange=[0,3]) :

	c = np.array([get_recovery_rate(hls['z_cl'][hls['snr_cl']>step], hls['mass'][hls['snr_cl']>step],
			speczRange, massRange, is_matched=(hls['mt_cross'][hls['snr_cl']>step] != None))['recovery'][0][0]
			for step in steps])
	cut = (cls['halo_mass']>massRange[0]) & (cls['sz']<=speczRange[1])
	p = np.array([get_recovery_rate(cls['z_cl'][cut & (cls['snr_cl']>step)], cls['mass'][cut & (cls['snr_cl']>step)],
			photozRange, richRange, is_matched=(cls['mt_cross'][cut & (cls['snr_cl']>step)] != None))['recovery'][0][0]
			for step in steps])
	t = np.array([sum((cls["snr_cl"]>step) & cut) / sum(cut) for step in steps])


	def find_nearest(array,value):
		idx = np.searchsorted(array, value, side="left")
		if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
			return idx-1
		else:
			return idx

	fig, axs = plt.subplots(1,1, figsize=figsize())
	axs.step(c, p, lw=1, c='k', ls='-', where='mid', label=f"$M>{massRange[0]:.2e}$, $z<{speczRange[1]:.1f}$ ({sum(cut)})")
	ind = [len(t) - find_nearest(t[::-1], perc) for perc in [0.01, 0.05, 0.1, 0.2, 0.4, 0.8, 0.99]]
	[axs.scatter(c[i], p[i], s=10, c='k', label=('(1, 5, 10, 20, 40, 80, 99)% of detections' if i==ind[0] else ''))
		for i in ind]
	axs.set_xlim(0,1.01)
	axs.set_ylim(0,1.01)
	axs.grid(which='major', lw=0.5)
	axs.grid(which='minor', lw=0.2)
	axs.legend(fontsize='small', frameon=True, loc='lower left', )
	axs.set_xlabel('Completeness', loc='right')
	axs.set_ylabel('Purity', loc='top')

	return fig, steps, c, p, t




#def completeness_skymap() :
#	info = recovery.skyplot(cat
