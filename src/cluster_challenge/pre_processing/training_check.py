import GCRCatalogs
from GCRCatalogs.helpers.base_filters import sample_filter, partition_filter
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.modeling import models, fitting
from astropy.table import Table
import numpy as np
import sys

#training file
table = ascii.read('/sps/lsst/groups/clusters/flexzboost/training_cosmoDC2/sample_cosmodc2_w10year_errors.dat')
table1 = table[table['i']<25]
table1 = table1[table1['i']>24]
table2 = table[table['i']<24]
table2 = table2[table2['i']>23]
table3 = table[table['i']<23]
table3 = table3[table3['i']>22]
table4 = table[table['i']<22]
table4 = table4[table4['i']>21]
table5 = table[table['i']<21]
#table5 = table5[table5['i']>20]
w2 = np.full(len(table2), len(table1)/len(table2))
w3 = np.full(len(table3), len(table1)/len(table3))
w4 = np.full(len(table4), len(table1)/len(table4))
w5 = np.full(len(table5), len(table1)/len(table5))
#print(table)

#flexzboost file
#file_1='/sps/lsst/users/tguillem/web/clusters/catalogs/cosmoDC2_photoz_flexzboost/v1/9559/galaxies.fits'
#file_1='/sps/lsst/users/tguillem/web/clusters/catalogs/cosmoDC2_photoz_flexzboost/v0/9559/galaxies.fits'
file_1='/sps/lsst/users/tguillem/web/clusters/catalogs/cosmoDC2_small_photoz_flexzboost/v2/9559/galaxies.fits'
galaxy_data_1 = Table.read(file_1)
#galaxy_data_1=galaxy_data_1[galaxy_data_1['mag_i_photoz']<25]
galaxy_data_2=galaxy_data_1[galaxy_data_1['mag_i']>26.5]
galaxy_data_3=galaxy_data_1[galaxy_data_1['mag_i']>25]
galaxy_data_3=galaxy_data_3[galaxy_data_3['mag_i']<26.5]
galaxy_data_4=galaxy_data_1[galaxy_data_1['mag_i']<25]
#galaxy_data_1=galaxy_data_1[galaxy_data_1['mag_i_photoz']<25]

file_2='/sps/lsst/users/tguillem/web/clusters/catalogs/cosmoDC2_photoz_flexzboost/v1/9554/galaxies.fits'
galaxy_data_5 = Table.read(file_2)
galaxy_data_5=galaxy_data_5[galaxy_data_5['photoz_mode']>0.01]
galaxy_data_5=galaxy_data_5[galaxy_data_5['photoz_mode']<2.99]
galaxy_data_5=galaxy_data_5[galaxy_data_5['mag_i']<25]
#galaxy_data_5=galaxy_data_5[galaxy_data_5['mag_i']>24]
filter_arr = []
for element in galaxy_data_5:
    if(abs(element['photoz_mode']-element['redshift'])<0.001):#*(1+element['redshift'])):
        filter_arr.append(True)
    else:
        filter_arr.append(False)
galaxy_data_5a = galaxy_data_5[filter_arr]
filter_arr = []
for element in galaxy_data_5:
    if(abs(element['photoz_mode']-element['redshift'])>0.001 and abs(element['photoz_mode']-element['redshift'])<0.01):#*(1+element['redshift'])):
        filter_arr.append(True)
    else:
        filter_arr.append(False)
galaxy_data_5b = galaxy_data_5[filter_arr]
filter_arr = []
for element in galaxy_data_5:
    if(abs(element['photoz_mode']-element['redshift'])>0.01):#*(1+element['redshift'])):
        filter_arr.append(True)
    else:
        filter_arr.append(False)
galaxy_data_5c = galaxy_data_5[filter_arr]

rescaling = len(table)/len(galaxy_data_1)
my_weights = np.full(len(galaxy_data_1), rescaling)
rescaling_2 = len(table)/len(galaxy_data_2)
my_weights_2 = np.full(len(galaxy_data_2), rescaling_2)
rescaling_3 = len(table)/len(galaxy_data_3)
my_weights_3 = np.full(len(galaxy_data_3), rescaling_3)
rescaling_4 = len(table)/len(galaxy_data_4)
my_weights_4 = np.full(len(galaxy_data_4), rescaling_4)
rescaling_5 = len(table)/len(galaxy_data_5)
my_weights_5 = np.full(len(galaxy_data_5), rescaling_5)

rescaling_5a = len(table)/len(galaxy_data_5a)
my_weights_5a = np.full(len(galaxy_data_5a), rescaling_5a)
rescaling_5b = len(table)/len(galaxy_data_5b)
my_weights_5b = np.full(len(galaxy_data_5b), rescaling_5b)
rescaling_5c = len(table)/len(galaxy_data_5c)
my_weights_5c = np.full(len(galaxy_data_5c), rescaling_5c)

nbins = 80
my_range = [0.8,1.6]
fig = plt.figure(figsize=(12,8))
#plt.hist(table1['sz'], density=False, range=my_range, bins = nbins, label="zspec 24-25", histtype='step', color = 'black');
#plt.hist(table2['sz'], density=False, range=my_range, bins = nbins, weights=w2, label="zspec 23-24", histtype='step', color = 'red');
plt.hist(table3['sz'], density=False, range=my_range, bins = nbins, weights=w3, label="zspec 22-23", histtype='step', color = 'blue');
#plt.hist(table4['sz'], density=False, range=my_range, bins = nbins, weights=w4, label="zspec 21-22", histtype='step', color = 'orange');
#plt.hist(table5['sz'], density=False, range=my_range, bins = nbins, weights=w5, label="zspec <21", histtype='step', color = 'purple');
#plt.hist(table['sz'], density=False, range=my_range, bins = nbins, label="zspec", histtype='step', color = 'black');
plt.hist(galaxy_data_5['photoz_mode'], density=False, range=my_range, bins = nbins, weights=my_weights_5, label="zmode", histtype='step', color = 'green');
#plt.hist(galaxy_data_5b['photoz_mode'], density=False, range=my_range, bins = nbins, weights=my_weights_5b, label="zmode 0.001<s<0.01", histtype='step', color = 'red');
#plt.hist(galaxy_data_5c['photoz_mode'], density=False, range=my_range, bins = nbins, weights=my_weights_5c, label="zmode s>0.01", histtype='step', color = 'blue');
#plt.hist(galaxy_data_1['redshift'], density=False, range=my_range, bins = nbins, weights=my_weights_2, label="ztrue hp1", histtype='step', color = 'blue');
#plt.hist(galaxy_data_1['redshift'], density=False, range=my_range, bins = nbins, weights=my_weights, label="zmode", histtype='step', color = 'red');
#plt.hist(data['photoz_mode'], density=False, range=my_range, bins = nbins, weights=my_weights, label="zmode", histtype='step', color = 'red');
#plt.hist(galaxy_data_2['redshift'], density=False, range=my_range, bins = nbins, weights=my_weights_2, label="ztrue hp1", histtype='step', color = 'blue');
#plt.hist(galaxy_data_3['redshift'], density=False, range=my_range, bins = nbins, weights=my_weights_3, label="ztrue hp1", histtype='step', color = 'red');
#plt.hist(galaxy_data_4['redshift'], density=False, range=my_range, bins = nbins, weights=my_weights_4, label="ztrue hp1", histtype='step', color = 'brown');
#plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
plt.xlabel("z");
plt.ylabel("galaxies");
#plt.xlim(1,1.5)
#plt.ylim(20,60)
#plt.legend([g1,g2], ["dr3", "cosmoDC2"], title = 'Photo-z', loc='upper right')
plt.legend(title = 'cosmoDC2 FlexZBoost', loc='upper right')
outpath = "/sps/lsst/users/tguillem/web/clusters/debug/flexzboost/"
plt.savefig(outpath+"redshift.png", bbox_inches='tight')
plt.close()

#mag correlation plot
plt.figure()
plt.scatter(galaxy_data_5a['mag_i'], galaxy_data_5a['mag_i_photoz'], marker='.',color = 'blue', s=10, alpha=0.3, label='galaxies')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([15, 27])
plt.ylim([15, 27])
plt.xlabel('mag_i')
plt.ylabel('mag_i_photoz')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'mag_mag.png', bbox_inches='tight')

#correlation plot
plt.figure()
plt.scatter(galaxy_data_1['redshift'], galaxy_data_1['mag_i_photoz'], marker='.',color = 'blue', s=0.1, alpha=0.3, label='galaxies')
#plt.xscale('log')
#plt.yscale('log')
plt.xlim([0, 3])
plt.ylim([15,30])
plt.xlabel('redshift')
plt.ylabel('mag_i_photoz')
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
plt.savefig(outpath+'mag_redshift.png', bbox_inches='tight')

#correlation plot
mags = ['mag_i','mag_r','mag_z','mag_y']
mags_photoz = ['mag_i_photoz','mag_r_photoz','mag_z_photoz','mag_y_photoz']
for i in range(4):
    plt.figure()
    plt.scatter(galaxy_data_1['redshift'], galaxy_data_1[mags[i]], marker='.',color = 'blue', s=0.1, alpha=0.3, label='galaxies')
    plt.xlim([0.3, 0.4])
    plt.ylim([16,30])
    plt.xlabel('redshift')
    plt.ylabel(mags[i])
    plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
    plt.savefig(outpath+mags[i]+'_redshift.png', bbox_inches='tight')
    
    plt.figure()
    plt.scatter(galaxy_data_1['redshift'], galaxy_data_1[mags_photoz[i]], marker='.',color = 'blue', s=1, alpha=0.3, label='galaxies')
    plt.xlim([0.3, 0.4])
    plt.ylim([16,30])
    plt.xlabel('redshift')
    plt.ylabel(mags_photoz[i])
    plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
    plt.savefig(outpath+mags_photoz[i]+'_redshift.png', bbox_inches='tight')

sys.exit()
