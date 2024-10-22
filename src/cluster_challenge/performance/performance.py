
import os, sys, yaml
import numpy as np

from clevar.catalog import ClCatalog
from clevar.match import get_matched_masks

import src.saving_figures as sfigs
from src.plotting_functions import make_bins, plot_hist
from src.pltFunc_matchingMetrics import make_ROC_curve, make_PurityCompleteness_curve



config = sys.argv[1]

# open config files
with open(config) as fstream:
	cfg = yaml.safe_load(fstream)

cat1 = cfg['cats']['cat1']
cat2 = cfg['cats']['cat2']
inpath = sfigs.make_path(cfg, base=cfg['paths']['performance']['in'])
outpath = sfigs.make_path(cfg, base=cfg['paths']['performance']['out'], addon='matching_metrics')


## read the catalogs.
cats = {cat1: ClCatalog.read_full(os.path.join(inpath, f"{cat1}.fits")),
	cat2: ClCatalog.read_full(os.path.join(inpath, f"{cat2}.fits"))}

z_bins    = make_bins('redshift', cfg, grain='both')
mass_bins = make_bins('mass', cfg, grain='both')
rich_bins = make_bins('richness', cfg, grain='both')

## share useful parameters between catalogs (SNR, halo_mass, richness, ...)
mt_mask_halo, mt_mask_cluster = get_matched_masks(cats[cat1], cats[cat2], 'cross')

cats[cat1]['snr_cl'] = 1e4
cats[cat1]['snr_cl'][mt_mask_halo] = cats[cat2]['snr'][mt_mask_cluster]

cats[cat1]['richness_cl'] = 1e6
cats[cat1]['richness_cl'][mt_mask_halo] = cats[cat2]['mass'][mt_mask_cluster]

cats[cat2]['halo_mass'] = 1e18
cats[cat2]['halo_mass'][mt_mask_cluster] = cats[cat1]['mass'][mt_mask_halo]

cats[cat2]['sz'] = np.zeros_like(cats[cat2]['mass'])
cats[cat2]['sz'][mt_mask_cluster] = cats[cat1]['z_cl'][mt_mask_halo]

cats[cat2]['halo_richness'] = 1e6
cats[cat2]['halo_richness'][mt_mask_cluster] = cats[cat1]['rich'][mt_mask_halo]

## make ROC curve
cut = (cats[cat2]['halo_mass'] > 1e14) & (cats[cat2]['sz']<1.5)
fig, TPR, FPR = make_ROC_curve(cats[cat2][cut])

saveas = 'ROCcurve'
sfigs.save_figure(fig, outpath, saveas=saveas)


## make Purity-Completeness curve
fig, snrSteps, c1, p1, t1 = make_PurityCompleteness_curve(cats[cat1], cats[cat2])
saveas = 'PurityCompleteness_Mgt1e14_zlt1.5'
sfigs.save_figure(fig, outpath, saveas=saveas)

fig, snrSteps, c2, p2, t2 = make_PurityCompleteness_curve(cats[cat1], cats[cat2], massRange=[10**13.5, 1e17])
saveas = 'PurityCompleteness_Mgt1e13.5_zlt1.5'
sfigs.save_figure(fig, outpath, saveas=saveas)

