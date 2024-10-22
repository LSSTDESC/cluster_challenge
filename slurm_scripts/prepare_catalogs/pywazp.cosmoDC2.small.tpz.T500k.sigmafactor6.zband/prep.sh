#!/bin/sh
#SBATCH --nodes=1
#SBATCH --job-name=prepare_catalogs
#SBATCH --time=2:00:00
#SBATCH --partition=lsst,htc
#SBATCH --ntasks=4
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/slurm_outputs/prepare_catalogs/pywazp.cosmoDC2.small.tpz.T500k.sigmafactor6.zband/prep.out
#SBATCH --mem=4GB
source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh
conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/src/cluster_challenge/prepare_catalogs/prepare_catalogs.py prepConfigs/pywazp.cosmoDC2.small.tpz.T500k.sigmafactor6.zband.cfg
