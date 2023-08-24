#!/bin/bash
#SBATCH --job-name=match_cluster_catalogs
#SBATCH --output=slurm_output.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=100:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=fail
##SBATCH --mail-type=end
##SBATCH --mail-user=solomon@lapp.in2p3.fr


source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh
conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/



## MATCH COSMODC2 AND WAZP.COSMODC2
python -u ../clevar_matching.py cosmoDC2_v1 wazp_cosmoDC2_v0 proximity
python -u ../clevar_matching.py cosmoDC2_v1 wazp_cosmoDC2_v1 proximity
python -u ../clevar_matching.py cosmoDC2_v1 wazp_cosmoDC2_v0 member
python -u ../clevar_matching.py cosmoDC2_v1 wazp_cosmoDC2_v1 member

## MATCH COSMODC2 AND WAZP.COSMODC2.FZB
python -u ../clevar_matching.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v0 proximity
python -u ../clevar_matching.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v1 proximity
python -u ../clevar_matching.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v0 member
python -u ../clevar_matching.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v1 member

## MATCH COSMODC2 AND WAZP.DC2
python -u ../clevar_matching.py cosmoDC2_v1 wazp_DC2_v0 proximity
python -u ../clevar_matching.py cosmoDC2_v1 wazp_DC2_v1 proximity

## MATCH COSMODC2 AND REDMAPPER.COSMODC2
python -u ../clevar_matching.py cosmoDC2_v1 redmapper_cosmoDC2_v0 proximity
python -u ../clevar_matching.py cosmoDC2_v1 redmapper_cosmoDC2_v1 proximity
python -u ../clevar_matching.py cosmoDC2_v1 redmapper_cosmoDC2_v0 member
python -u ../clevar_matching.py cosmoDC2_v1 redmapper_cosmoDC2_v1 member

## MATCH COSMODC2 AND REDMAPPER.DC2
python -u ../clevar_matching.py cosmoDC2_v1 redmapper_DC2_v0 proximity
python -u ../clevar_matching.py cosmoDC2_v1 redmapper_DC2_v1 proximity

## MATCH WAZP.COSMODC2 AND WAZP.COSMODC2.FZB
python -u ../clevar_matching.py wazp_cosmoDC2_v0 wazp_cosmoDC2.fzb_v0 proximity
python -u ../clevar_matching.py wazp_cosmoDC2_v1 wazp_cosmoDC2.fzb_v1 proximity
python -u ../clevar_matching.py wazp_cosmoDC2_v0 wazp_cosmoDC2.fzb_v0 member
python -u ../clevar_matching.py wazp_cosmoDC2_v1 wazp_cosmoDC2.fzb_v1 member

## MATCH WAZP.COSMODC2 AND WAZP.DC2
python -u ../clevar_matching.py wazp_cosmoDC2_v0 wazp_DC2_v0 proximity
python -u ../clevar_matching.py wazp_cosmoDC2_v1 wazp_DC2_v1 proximity

## MATCH WAZP.COSMODC2.FZB AND WAZP.DC2
python -u ../clevar_matching.py wazp_cosmoDC2.fzb_v0 wazp_DC2_v0 proximity
python -u ../clevar_matching.py wazp_cosmoDC2.fzb_v1 wazp_DC2_v1 proximity

## MATCH REDMAPPER.COSMODC2 AND WAZP.COSMODC2
python -u ../clevar_matching.py redmapper_cosmoDC2_v0 wazp_cosmoDC2_v0 proximity
python -u ../clevar_matching.py redmapper_cosmoDC2_v1 wazp_cosmoDC2_v1 proximity
python -u ../clevar_matching.py redmapper_cosmoDC2_v0 wazp_cosmoDC2_v0 member
python -u ../clevar_matching.py redmapper_cosmoDC2_v1 wazp_cosmoDC2_v1 member

## MATCH REDMAPPER.DC2 AND WAZP.DC2
python -u ../clevar_matching.py redmapper_DC2_v0 wazp_DC2_v0 proximity
python -u ../clevar_matching.py redmapper_DC2_v1 wazp_DC2_v1 proximity

### MATCH COSMODC2 AND AMICO.COSMODC2.FZB
#python -u ../clevar_matching.py cosmoDC2_v1 amico_cosmoDC2.fzb_v0 proximity
#python -u ../clevar_matching.py cosmoDC2_v1 amico_cosmoDC2.fzb_v0 member
