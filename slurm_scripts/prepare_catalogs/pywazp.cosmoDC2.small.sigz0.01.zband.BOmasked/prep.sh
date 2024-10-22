#!/bin/sh
#SBATCH --nodes=1
#SBATCH --job-name=prepare_catalogs
#SBATCH --time=2:00:00
#SBATCH --partition=lsst,htc
#SBATCH --ntasks=4
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/slurm_outputs/prepare_catalogs/pywazp.cosmoDC2.small.sigz0.01.zband.BOmasked/prep.out
#SBATCH --mem=4GB
module load conda
conda activate /sps/lsst/users/rsolomon/conda_envs/desc_dev/
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/src/cluster_challenge/prepare_catalogs/prepare_catalogs.py prepConfigs/pywazp.cosmoDC2.small.sigz0.01.zband.BOmasked.cfg
