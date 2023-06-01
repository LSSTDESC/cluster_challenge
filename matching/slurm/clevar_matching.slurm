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
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py wazp cosmoDC2 v0 member
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py wazp cosmoDC2 v1 member
## match cosmoDC2 and wazp.DC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py wazp DC2 v0 member


## match cosmoDC2 and redmapper.cosmoDC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py redmapper cosmoDC2 v0 member
## match cosmoDC2 and redmapper.DC2
python -u /pbs/throng/lsst/users/rsolomon/cluster_challenge/matching/clevar_matching.py redmapper DC2 v0 member
