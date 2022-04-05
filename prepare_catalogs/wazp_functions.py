import numpy as np
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
from astropy.table import Table
from astropy.io import ascii
import pickle

#######################################################
#function to open, read and select info in WaZP catalog
#######################################################

def wazp_cat_open(wazp_cat_name, wazp_members_cat_name, min_richness):
    cat_wazp = Table.read(wazp_cat_name)
    cat_wazp_members = Table.read(wazp_members_cat_name)
    #print(cat_wazp)
    
    #cat_wazp
    cat_wazp.rename_column('RA', 'ra')
    cat_wazp.rename_column('DEC', 'dec')
    cat_wazp.rename_column('zp_bright', 'redshift')
    ###DES
    #cat_wazp.rename_column('ZP_bright', 'redshift')
    
    #cat_wazp_members
    cat_wazp_members.rename_column('RA', 'ra')
    cat_wazp_members.rename_column('DEC', 'dec')
    cat_wazp_members.rename_column('ZP', 'redshift')
    #cat_wazp_members.rename_column('ZP_CL', 'redshift')

    #filter for richness
    max_richness = 1000
    filter = np.logical_and(cat_wazp['NGALS']>=min_richness,cat_wazp['NGALS']<max_richness)
    #filter_2 = np.logical_and(cat_wazp['redshift']>=0.5,cat_wazp['redshift']<0.6) 
    #filter = np.logical_and(filter_1,filter_2)
    cat_wazp_filter = cat_wazp[filter]
    filter =  np.logical_and(cat_wazp_members['NGALS']>min_richness,cat_wazp_members['NGALS']<max_richness)
    cat_wazp_members_filter = cat_wazp_members[filter]
    return cat_wazp_filter, cat_wazp_members_filter
    
    #return cat_wazp, cat_wazp_members 

def mstar_i(redshift):
    bin = round(redshift-0.020,2)/0.010
    bin_r = round(bin) +1
    if(bin_r>248):
         bin_r = 248
    return bin_r

def mstar_z(redshift):
    bin = round(redshift,2)/0.05
    bin_r = round(bin)
    if(bin_r>55):
         bin_r = 55
    return bin_r

def clean_2D(input_file):

    my2D = pickle.load(open(input_file,"rb")) 
    
    h = np.array(my2D[0])
    xedges = np.array(my2D[1])
    yedges = np.array(my2D[2])
    
    nbins_x = xedges.size-1
    nbins_y = yedges.size-1

    #set nan to 0
    for ix in range(nbins_x):
        for iy in range(nbins_y):
            if np.isnan(h[ix,iy]):
                h[ix,iy]=0

    return h, xedges, yedges, nbins_x, nbins_y
