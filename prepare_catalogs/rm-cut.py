import numpy as np
import matplotlib.pyplot as plt
from clevar.match import output_matched_catalog
from clevar.catalog import ClCatalog

mass_min = 17

rm_path = '/sps/lsst/users/tguillem/DESC/desc_april_2022/cluster_challenge/clevar_catalogs/redmapper/full_pmem/cosmoDC2_v1.1.4_redmapper_v0.8.1/'

#Load cat
rm = ClCatalog.read(rm_path + 'Catalog.fits', 'rm', full =True)
rm.read_members(rm_path+'Catalog_members.fits', full =True)

#Cut cat
rm = rm[rm['mass']>=mass_min]

out = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/RedMapper/'
rm.write(out + 'Catalog_ls17' + '.fits', overwrite = True)
rm.members.write(out + 'Catalog_members_ls17' + '.fits', overwrite = True)