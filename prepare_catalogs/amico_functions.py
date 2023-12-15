import numpy as np
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
from astropy.table import Table
from astropy.io import ascii
import pickle

#######################################################
# function to open, read and select info in AMICO catalog
# ######################################################

def mstar_i(redshift):
    bin = round(redshift,2)/0.010
    bin_r = round(bin) +1
    if(bin_r>55):
         bin_r = 55
    return bin_r

def mstar_z(redshift):
    bin = round(redshift,2)/0.05
    bin_r = round(bin)
    if(bin_r>55):
         bin_r = 55
    return bin_r

def mstar_y(redshift):
    bin = round(redshift,2)/0.010
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
