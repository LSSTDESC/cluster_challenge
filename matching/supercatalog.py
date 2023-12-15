import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
from clevar.catalog import ClCatalog
from numpy.ma import masked
import sys
from clevar.match import output_matched_catalog

p_matching = True
mb_matching = False

if p_matching is True:
    matching = 'p'
elif mb_matching is True:
    matching = 'mb'  
else :
    print('invalid matching')


inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/'
#am-cdc
matching_path = inpath + 'amico_cosmoDC2/mag_y/'
cam1 = ClCatalog.read_full(matching_path + 'c1_' + matching + '_ls4512z115.fits')
cdc1 = ClCatalog.read_full(matching_path + 'c2_'+ matching +'_ls4512z115.fits')
print('am & cdc loaded')
output_matched_catalog(matching_path + 'c1_' + matching + '_ls4512z115.fits', matching_path + 'c2_' + matching + '_ls4512z115.fits'
                       ,matching_path + 'output_catalog_' + matching +'_ls4512z115.fits', cam1, cdc1, matching_type='cross', overwrite=True)
c_merged_13 = Table.read(matching_path+'output_catalog_' + matching + '_ls4512z115.fits')
print('c_merged_13 ok')


#am-rm
matching_path = inpath + 'amico_redmapper/mag_y/'
cam2 = ClCatalog.read_full(matching_path + 'c1_' + matching + '_ls4512z115.fits')
crm1 = ClCatalog.read_full(matching_path + 'c2_' + matching + '_ls4512z115.fits')
print('am & rm loaded')
output_matched_catalog(matching_path + 'c1_' + matching + '_ls4512z115.fits', matching_path + 'c2_' + matching + '_ls4512z115.fits'
                       ,matching_path+'output_catalog_' + matching + '_ls4512z115.fits', cam2, crm1, matching_type='cross', overwrite=True)
c_merged_12 = Table.read(matching_path+'output_catalog_' + matching + '_ls4512z115.fits')
print('c_merged_12 ok')


#rm-cdc2
matching_path = inpath + 'redmapper_cosmoDC2/'
crm2 = ClCatalog.read_full(matching_path + 'c1_' + matching + '_ls4512z115.fits')
cdc2 = ClCatalog.read_full(matching_path + 'c2_' + matching + '_ls4512z115.fits')
print('rm & cdc loaded')
output_matched_catalog(matching_path + 'c1_' + matching + '_ls4512z115.fits', matching_path+'c2_' + matching + '_ls4512z115.fits'
                       ,matching_path+'output_catalog_' + matching +'_ls4512z115.fits', crm2, cdc2, matching_type='cross', overwrite=True)
c_merged_23 = Table.read(matching_path+'output_catalog_' + matching + '_ls4512z115.fits')
print('c_merged_23 ok')

#Merging all

def rename(cat, c1name1, c1rename1, c1name2, c1rename2, match): 
    names = cat.colnames
    rename = []
    print(names)
    for i in range(len(names)):
        if c1name1 in names[i]:
            rename.append(c1rename1 + names[i][len(c1name1):])
        elif c1name2 in names[i]: 
            rename.append(c1rename2 + names[i][len(c1name2):])
        elif names[i] == 'id':
            rename.append('id_' + match)
    print(rename)
    return rename

l = [len(c_merged_12), len(c_merged_13), len(c_merged_23)]
r = max(l)

def shaping(cat, renamed): #renamed = rename(cat, c1name1, c1rename1, c1name2, c1rename2, match) see example bellow
    if len(cat)!=r:
        rename_l = renamed
        cat.rename_columns(cat.colnames, rename_l)
        Z_arr = np.zeros([r-len(cat),len(rename_l)])
        Z = Table(Z_arr, names = rename_l)
        for name in rename_l:
            if 'mt' in name:
                Z[name] = Z[name].astype(bytes)
                cat[name] = cat[name].astype(bytes)
            elif 'id' in name:
                Z[name] = Z[name].astype(bytes)
                cat[name] = cat[name].astype(bytes)
        stacked = vstack([cat, Z])
        return stacked
                         
    else :
        rename_l = renamed
        cat.rename_columns(cat.colnames, rename_l)
        return cat
print('shaping 12')    
c_merged_12_renamed = rename(c_merged_12, 'cat1', 'cat12-1', 'cat2', 'cat12-2', '12')
c_merged_12 = shaping(c_merged_12, c_merged_12_renamed)
print('shaping 13') 
c_merged_13_renamed = rename(c_merged_13, 'cat1', 'cat13-1', 'cat2', 'cat13-3', '13')
c_merged_13 = shaping(c_merged_13, c_merged_13_renamed)
print('shaping 23') 
c_merged_23_renamed = rename(c_merged_23, 'cat1', 'cat23-2', 'cat2', 'cat23-3', '23')
c_merged_23 = shaping(c_merged_23, c_merged_23_renamed)

supercatalog = hstack([c_merged_12,c_merged_13])
supercatalog = hstack([supercatalog,c_merged_23])
supercatalog.write('/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/' + matching + 'supercatalog_ls4512z115_my.fits', overwrite = True)