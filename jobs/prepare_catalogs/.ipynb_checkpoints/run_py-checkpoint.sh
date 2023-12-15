#!/bin/sh
echo "DESC setup"
source /pbs/home/n/namourou/test_jupyter/cluster_challenge/jobs/setup.sh
#conda env list

echo "GO1"

###CATALOGS

python /pbs/home/n/namourou/test_jupyter/cluster_challenge/prepare_catalogs/merge_amico_mb.py

echo "G02"

python /pbs/home/n/namourou/test_jupyter/cluster_challenge/prepare_catalogs/read_amico.py

echo "GO3"
