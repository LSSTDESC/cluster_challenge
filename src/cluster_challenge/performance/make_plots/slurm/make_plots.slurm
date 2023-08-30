#!/bin/bash
#SBATCH --job-name=plot_matched_cluster_catalogs
##SBATCH --output=slurm_output.out
#SBATCH --output=/dev/null --error=/dev/null
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --time=10:00:00
##SBATCH --mail-type=begin
##SBATCH --mail-type=fail
##SBATCH --mail-type=end
##SBATCH --mail-user=solomon@lapp.in2p3.fr


source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh
conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/

if [ -n $SLURM_JOB_ID ]; then
        SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
else
        SCRIPT_PATH=$(realpath $0)
fi

SHARED_PATH=$(dirname  $SCRIPT_PATH)

exec > $SHARED_PATH/slurm_output.out

echo $(date -u +%Y%m%d-%H:%M:%S)


## PLOT COSMODC2 AND WAZP.COSMODC2
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 wazp_cosmoDC2_v0 proximity
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 wazp_cosmoDC2_v1 proximity
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 wazp_cosmoDC2_v0 member
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 wazp_cosmoDC2_v1 member

## PLOT COSMODC2 AND WAZP.COSMODC2.FZB
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v0 proximity
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v1 proximity
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v0 member
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 wazp_cosmoDC2.fzb_v1 member

## PLOT COSMODC2 AND WAZP.DC2
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 wazp_DC2_v0 proximity
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 wazp_DC2_v1 proximity

## PLOT COSMODC2 AND REDMAPPER.COSMODC2
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 redmapper_cosmoDC2_v0 proximity
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 redmapper_cosmoDC2_v1 proximity
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 redmapper_cosmoDC2_v0 member
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 redmapper_cosmoDC2_v1 member

## PLOT COSMODC2 AND REDMAPPER.DC2
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 redmapper_DC2_v0 proximity
python -u $SHARED_PATH/../make_plots.py cosmoDC2_v1 redmapper_DC2_v1 proximity

## PLOT WAZP.COSMODC2 AND WAZP.COSMODC2.FZB
python -u $SHARED_PATH/../make_plots.py wazp_cosmoDC2_v0 wazp_cosmoDC2.fzb_v0 proximity
python -u $SHARED_PATH/../make_plots.py wazp_cosmoDC2_v1 wazp_cosmoDC2.fzb_v1 proximity
python -u $SHARED_PATH/../make_plots.py wazp_cosmoDC2_v0 wazp_cosmoDC2.fzb_v0 member
python -u $SHARED_PATH/../make_plots.py wazp_cosmoDC2_v1 wazp_cosmoDC2.fzb_v1 member

## PLOT WAZP.COSMODC2 AND WAZP.DC2
python -u $SHARED_PATH/../make_plots.py wazp_cosmoDC2_v0 wazp_DC2_v0 proximity
python -u $SHARED_PATH/../make_plots.py wazp_cosmoDC2_v1 wazp_DC2_v1 proximity

## PLOT WAZP.COSMODC2.FZB AND WAZP.DC2
python -u $SHARED_PATH/../make_plots.py wazp_cosmoDC2.fzb_v0 wazp_DC2_v0 proximity
python -u $SHARED_PATH/../make_plots.py wazp_cosmoDC2.fzb_v1 wazp_DC2_v1 proximity

## PLOT REDMAPPER.COSMODC2 AND WAZP.COSMODC2
python -u $SHARED_PATH/../make_plots.py redmapper_cosmoDC2_v0 wazp_cosmoDC2_v0 proximity
python -u $SHARED_PATH/../make_plots.py redmapper_cosmoDC2_v1 wazp_cosmoDC2_v1 proximity
python -u $SHARED_PATH/../make_plots.py redmapper_cosmoDC2_v0 wazp_cosmoDC2_v0 member
python -u $SHARED_PATH/../make_plots.py redmapper_cosmoDC2_v1 wazp_cosmoDC2_v1 member

## PLOT REDMAPPER.DC2 AND WAZP.DC2
python -u $SHARED_PATH/../make_plots.py redmapper_DC2_v0 wazp_DC2_v0 proximity
python -u $SHARED_PATH/../make_plots.py redmapper_DC2_v1 wazp_DC2_v1 proximity


echo $(date -u +%Y%m%d-%H:%M:%S)

