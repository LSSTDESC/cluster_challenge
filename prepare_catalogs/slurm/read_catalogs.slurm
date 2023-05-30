#!/bin/bash
#SBATCH --job-name=reduce_cluster_catalogs
#SBATCH --output=/pbs/throng/lsst/users/rsolomon/cluster_challenge/prepare_catalogs/slurm/slurm_output.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=100:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=fail
#SBATCH --mail-type=end
#SBATCH --mail-user=solomon@lapp.in2p3.fr


. /pbs/home/r/rsolomon/.s

#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/prepare_catalogs/read_cosmoDC2.py v0
#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/prepare_catalogs/read_cosmoDC2.py v1

#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/prepare_catalogs/read_wazp.cosmoDC2.py v0
#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/prepare_catalogs/read_wazp.cosmoDC2.py v1

#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/prepare_catalogs/read_wazp.DC2.py v0

#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/prepare_catalogs/read_redmapper.cosmoDC2.py v0

#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/prepare_catalogs/read_redmapper.DC2.py v0

#python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/prepare_catalogs/read_amico.cosmoDC2.py v0
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/prepare_catalogs/read_amico.cosmoDC2.py v1
