#!/usr/bin/env python
# coding: utf-8
print('Merging amico files')
###import
import numpy as np
from astropy.table import Table, vstack
import sys
import os


# ---Paths---#
###CosmoDC2
####magi
inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/'
####magy
#inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_y/'

###DC2
####magy
#inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/DC2/mag_y/'

outpath = inpath

print('Reading files at', inpath, 'which will be saved at', outpath)


# ---Setting tiles up & reading/merging tables---#

#CosmoDC2
tiles = ['10066', '10067', '10068', '10069', '10070', '10071', '10072', '10073', '10074', '10193', '10194', '10195', '10196', '10197', '10198', '10199', '10200', '10201', '10202', '10321', '10322', '10323', '10324', '10325', '10326', '10327', '10328', '10329', '10444', '10445', '10446', '10447', '10448', '10449', '10450', '10451', '10452', '8786', '8787', '8788', '8789', '8790', '8791', '8792', '8793', '8794', '8913', '8914', '8915', '8916', '8917', '8918', '8919', '8920', '8921', '9042', '9043', '9044', '9045', '9046', '9047', '9048', '9049', '9050', '9169', '9170', '9171', '9172', '9173', '9174', '9175', '9176', '9177', '9178', '9298', '9299', '9300', '9301', '9302', '9303', '9304', '9305', '9306', '9425', '9426', '9427', '9428', '9429', '9430', '9431', '9432', '9433', '9434', '9554', '9555', '9556', '9557', '9558', '9559', '9560', '9561', '9562', '9681', '9682', '9683', '9684', '9685', '9686', '9687', '9688', '9689', '9690', '9810', '9811', '9812', '9813', '9814', '9815', '9816', '9817', '9818', '9937', '9938', '9939', '9940', '9941', '9942', '9943', '9944', '9945', '9946']

#DC2
#tiles = [5074, 5073, 5072, 5071, 5070, 5069, 5068, 5067, 5066, 5065, 4860, 4859, 4858, 4857, 4856, 4855, 4854, 4853, 4852, 4851, 4850, 4648, 4647, 4646, 4645, 4644, 4643, 4642, 4641, 4640, 4639, 4638, 4637, 4636, 4441, 4440, 4439, 4438, 4437, 4436, 4435, 4434, 4433, 4432, 4431, 4430, 4429, 4236, 4235, 4234, 4233, 4232, 4231, 4230, 4229, 4228, 4227, 4226, 4225, 4224, 4035, 4034, 4033, 4032, 4031, 4030, 4029, 4028, 4027, 4026, 4025, 4024, 4023, 3837, 3836, 3835, 3834, 3833, 3832, 3831, 3830, 3829, 3828, 3827, 3826, 3825, 3643, 3642, 3641, 3640, 3639, 3638, 3637, 3636, 3635, 3634, 3633, 3632, 3631, 3453, 3452, 3451, 3450, 3449, 3448, 3447, 3446, 3445, 3444, 3443, 3442, 3441, 3268, 3267, 3266, 3265, 3264, 3263, 3262, 3261, 3260, 3259, 3258, 3257, 3256, 3086, 3085, 3084, 3083, 3082, 3081, 3080, 3079, 3078, 3077, 3076, 3075, 3074, 2908, 2907, 2906, 2905, 2904, 2903, 2902, 2901, 2900, 2899, 2898, 2897, 2896]

large_amico_data = Table.read(inpath + str(tiles[0]) + '_map_associations.fits') 

for tile in tiles[1:] : 
    print('processing tile =', tile)
    amico_data = Table.read(inpath + str(tile) + '_map_associations.fits')
    large_amico_data = vstack([large_amico_data, amico_data])
    print('Done')


# ---Saving---#

large_amico_data.write(outpath + 'all_maps.fits', overwrite=True)

sys.exit()