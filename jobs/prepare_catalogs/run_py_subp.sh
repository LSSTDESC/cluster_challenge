#!/bin/sh
echo "DESC setup"
source /pbs/home/n/namourou/test_jupyter/cluster_challenge/jobs/setup.sh
#conda env list

echo "GO1"

###CATALOGS
python /pbs/home/n/namourou/test_jupyter/cluster_challenge/prepare_catalogs/amico_add_mag.py ${tile}

echo "G02"

python /pbs/home/n/namourou/test_jupyter/cluster_challenge/prepare_catalogs/read_amico_mb_subp.py ${tile}

echo "GO3"

