#!/bin/bash
#SBATCH --job-name=plot_matched_cluster_catalogs
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/slurm/slurm_output.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --time=10:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=solomon@lapp.in2p3.fr


. /pbs/home/r/rsolomon/.s



## PLOT COSMODC2 AND WAZP.COSMODC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 wazp_cosmoDC2_v0 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 wazp_cosmoDC2_v1 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 wazp_cosmoDC2_v0 member
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 wazp_cosmoDC2_v1 member

## PLOT COSMODC2 AND WAZP.COSMODC2.FZB
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v0 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v1 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v0 member
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v1 member

## PLOT COSMODC2 AND WAZP.DC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 wazp_DC2_v0 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 wazp_DC2_v1 proximity

## PLOT COSMODC2 AND REDMAPPER.COSMODC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 redmapper_cosmoDC2_v0 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 redmapper_cosmoDC2_v1 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 redmapper_cosmoDC2_v0 member
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 redmapper_cosmoDC2_v1 member

## PLOT COSMODC2 AND REDMAPPER.DC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 redmapper_DC2_v0 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py cosmoDC2_v1 redmapper_DC2_v1 proximity

## PLOT WAZP.COSMODC2 AND WAZP.DC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py wazp_cosmoDC2_v0 wazp_DC2_v0 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py wazp_cosmoDC2_v1 wazp_DC2_v1 proximity

## PLOT WAZP.COSMODC2.FZB AND WAZP.DC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py wazp_cosmoDC2.fzb_v0 wazp_DC2_v0 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py wazp_cosmoDC2.fzb_v1 wazp_DC2_v1 proximity

## PLOT REDMAPPER.COSMODC2 AND WAZP.COSMODC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py redmapper_cosmoDC2_v0 wazp_cosmoDC2_v0 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py redmapper_cosmoDC2_v1 wazp_cosmoDC2_v1 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py redmapper_cosmoDC2_v0 wazp_cosmoDC2_v0 member
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py redmapper_cosmoDC2_v1 wazp_cosmoDC2_v1 member

## PLOT REDMAPPER.DC2 AND WAZP.DC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py redmapper_DC2_v0 wazp_DC2_v0 proximity
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/performance/make_plots/make_plots.py redmapper_DC2_v1 wazp_DC2_v1 proximity

