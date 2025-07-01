
from mt_utils import mtCatalogs
import sys


config = sys.argv[1]

mtcats = mtCatalogs(config)


print(f'clevar version: {mtcats.get_clevar_version()}')

print("getting catalogs", end='...')
mtcats.get_catalogs()
print('DONE')

print("reading members", end='...')
mtcats.get_mb_catalogs()
print('DONE')

print('applying footprints', end='...')
mtcats.apply_mt_footprints_mask()
print('DONE')

print('perform matching', end='...')
mtcats.match_catalogs()
print('DONE')

print('writing matched files', end='...')
mtcats.write_mtCatalogs()
print('DONE')
