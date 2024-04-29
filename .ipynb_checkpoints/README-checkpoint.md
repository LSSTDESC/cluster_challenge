# cluster_challenge
Comparison of cluster algorithm performances on cosmoDC2 and DC2 datasets   
Developed for the DC2 project: https://portal.lsstdesc.org/DESCPub/app/PB/show_project?pid=248   
Common framework for the analysis of cosmoDC2, redMaPPer, WaZP and AMICO catalogs   
Using CLEVAR package   
For CLEVAR installation: pip install -e clevar  (after python setup.py install --user)

## DESC environment
### Create a DESC conda environment (to be done only once)   
source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh    
conda create --clone desc -p desc_july_2022_v0   
conda activate desc_july_2022_v0   
conda install mysql-connector-python -c conda-forge   
### Setup (to be done every time you log in)  
source /pbs/throng/lsst/software/desc/common/miniconda/setup_current_python.sh   
conda activate desc_july_2022_v0   
