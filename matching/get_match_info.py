import numpy as np
from astropy.table import Table
from clevar.catalog import ClCatalog
from numpy.ma import masked
import sys

#---Parameters selection---#

inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/'
cat_path = inpath + str(sys.argv[1])

print('Using catalogs located at :', cat_path)

##---Select catalogs to match---##

outpath = cat_path

print('Saving plots at :', outpath)

##---Select matching---##
if str(sys.argv[2]) == 'p_matching':
    matching = 'p'
elif str(sys.argv[2]) == 'mb_matching':
    matching = 'mb'
    
#---Reading Catlogs---#
c1 = Table.read(cat_path + 'c1_' + matching + '.fits')
c2 = Table.read(cat_path + 'c2_' + matching + '.fits')

#---Defining function to read matching characteristics---#
def get_info(cat):
    j, k, l, m, n = 0, 0, 0, 0, 0
    for i in range(len(cat['mt_self'])):
        if cat['mt_self'][i] is not masked:
            j+=1
        if cat['mt_other'][i] is not masked:
            k+=1
        if cat['mt_multi_self'][i] is not masked:
            l+=1
        if cat['mt_multi_other'][i] is not masked:
            m+=1
        if cat['mt_cross'][i] is not masked:
            n+=1  
    return j, k, l, m, n

##---Use function and print information---##
j, k, l, m, n = get_info(c1)
print('c1 infos from ' + matching + '-matching with c2','\n#c1', '\ntotal objects', len(c1), '\nunique (self):', j, '\nunique (other):', k, '\nmultiple (self):', l, '\nmultiple (other):', m, '\ncross:', n)
j2, k2, l2, m2, n2 = get_info(c2)
print('\nc2 infos from ' + matching + '-matching with c1','\n#c2', '\ntotal objects', len(c2), '\nunique (self):', j2, '\nunique (other):', k2, '\nmultiple (self):', l2, '\nmultiple (other):', m2, '\ncross:', n2)

###---Reunite informations in a table & saves it---###
dic = {'id' : ['c1', 'c2'] , 'total' : [len(c1),len(c2)] , 'unique (self)' : [j, j2], 'unique (other)' : [k, k2], 'multiple (self)' : [l, l2], 'multiple (other)' : [m, m2], 'cross' : [n, n2]}
t = Table(dic)
print(t)
t.write(outpath + 'match_info_' + matching + '.fits', overwrite = True)