#!/bin/sh
echo "DESC setup"
source /pbs/home/n/namourou/test_jupyter/cluster_challenge/jobs/setup.sh
#conda env list

echo "GO1"
#outpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/" 

# ##CATALOGS
echo 'charged:' ${param1} ${param2} ${param3}
python /pbs/home/n/namourou/test_jupyter/cluster_challenge/matching/clevar_matching.py ${param1} ${param2} ${param3}

echo "G02"

#python /pbs/home/n/namourou/test_jupyter/cluster_challenge/matching/get_match_info.py /sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/${param}/mag_i/

echo "GO3"
