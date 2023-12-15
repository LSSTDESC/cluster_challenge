#!/bin/sh
echo "DESC setup"
source /pbs/home/n/namourou/test_jupyter/cluster_challenge/jobs/setup.sh
#conda env list


#outpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/" 

# ##CATALOGS
#echo 'charged:' ${param1} ${param2}
#python /pbs/home/n/namourou/test_jupyter/cluster_challenge/notebooks/NB_analysis/am-cdc/magi-magr.py
echo "matching infos"
#python /pbs/home/n/namourou/test_jupyter/cluster_challenge/matching/get_match_info.py ${param1} ${param2}

echo "completeness"
#python /pbs/home/n/namourou/test_jupyter/cluster_challenge/performance/completeness_m.py ${param1} ${param2}
#python /pbs/home/n/namourou/test_jupyter/cluster_challenge/performance/completeness_z.py ${param1} ${param2}
#python /pbs/home/n/namourou/test_jupyter/cluster_challenge/performance/completeness_2D.py ${param1} ${param2}

echo "purity"

#python /pbs/home/n/namourou/test_jupyter/cluster_challenge/performance/purity.py ${param1} ${param2}

echo "overmerging & fragmentation"

python /pbs/home/n/namourou/test_jupyter/cluster_challenge/performance/overmerging.py ${param1} ${param2}
#python /pbs/home/n/namourou/test_jupyter/cluster_challenge/performance/fragmentation.py ${param1} ${param2}

