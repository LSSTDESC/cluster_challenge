#!/bin/bash
#SBATCH --job-name=match_cluster_catalogs
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/slurm/slurm_output.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=100:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=solomon@lapp.in2p3.fr


. /pbs/home/r/rsolomon/.s


## match cosmoDC2 and wazp.cosmoDC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py cosmoDC2 v1 wazp_cosmoDC2 v0 member_matching
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py cosmoDC2 v1 wazp_cosmoDC2 v1 member_matching

## match cosmoDC2 and redmapper.cosmoDC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py cosmoDC2 v1 redmapper_cosmoDC2 v0 member_matching

## match redmapper.cosmoDC2 and redmapper.DC2
#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py redmapper_cosmoDC2 v0 redmapper_DC2 v0 member_matching

## match wazp.cosmoDC2 and wazp.DC2
#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py wazp_cosmoDC2 v0 wazp_DC2 v0 member_matching
#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py wazp_cosmoDC2 v1 wazp_DC2 v0 member_matching

## match redmapper.cosmoDC2 and wazp.cosmoDC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py redmapper_cosmoDC2 v0 wazp_cosmoDC2 v0 member_matching
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py redmapper_cosmoDC2 v0 wazp_cosmoDC2 v1 member_matching

## match redmapper.DC2 and wazp.DC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py wazp_DC2 v0 redmapper_DC2 v0 member_matching
