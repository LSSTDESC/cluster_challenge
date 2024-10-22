#!/bin/sh
#SBATCH --nodes=1
#SBATCH --job-name=prepare_catalogs
#SBATCH --time=2:00:00
#SBATCH --partition=lsst,htc
#SBATCH --ntasks=4
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/slurm_outputs/prepare_catalogs/pywazp.DC2.tpz.T500k.zband/prep.out
#SBATCH --mem=8GB
module load conda
conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/src/cluster_challenge/prepare_catalogs/prepare_catalogs.py prepConfigs/pywazp.DC2.small.tpz.T500k.zband.cfg
