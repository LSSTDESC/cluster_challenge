import numpy as np
import glob
import os
import pickle
import struct

import matplotlib.pyplot as plt




def make_path(cat1, cat2, mt_method, mt_pref, mt_params, base, addon=None) :

	algo1 = '.'.join(cat1.split('.')[:-1])
	algo2 = '.'.join(cat2.split('.')[:-1])

	v1 = cat1.split('.')[-1]	
	v2 = cat2.split('.')[-1]	
	
	path  = base
	path += f"{algo1}_{algo2}/"	## EXAMPLE CATS: cosmoDC2.v0 wazp.cosmoDC2.fzb.v0
	path += f"{v1}_{v2}/"		## RESULTING PATH: /cosmoDC2_wazp.cosmoDC2.fzb/v0_v0/
	
	if mt_method == 'proximity' :
		path += f"proximity_matching/deltaz_{mt_params[0]}_matchradius_{mt_params[1]}mpc"
	else :
		path += f"member_matching/fshare_{mt_params[0]}"

	path += f"_pref_{mt_pref}/"
	
	if addon is None :
		return path
	else :
		path += f"{addon}/"
		return path



def save_figure(fig, outpath, saveas, png=True, pdf=False, pkl=False) :
	ftypes = ['png', 'pdf', '.pickle']
	outpaths = np.array([f"{outpath}{ftype}/" for ftype in ftypes])[[png, pdf, pkl]]
	for path in outpaths :
		if not os.path.exists(path) :
			print(f"Making directory:\t{path}")
			os.makedirs(path)
	if png :
		fig.savefig(outpaths[0] + saveas + f".{ftypes[0]}", bbox_inches='tight')
	if pdf :
		fig.savefig(outpaths[1] + saveas + f".{ftypes[1]}", bbox_inches='tight')
	if pkl :
		pickle.dump(fig, open(outpaths[2] + saveas + f"{ftypes[2]}", 'wb'))
	plt.close(fig)
	
	update_index_file_for_html_display(outpath)



def update_index_file_for_html_display(outpath, description='') :
	files = list(filter(os.path.isfile, glob.glob(f'{outpath}png/*')))
	files.sort(key=lambda x: os.path.getmtime(x))
	
	width = []
	height = []
	for f in files :
		with open(f, 'rb') as fhandle :
			head = fhandle.read(24)
			w, h = struct.unpack('>ii', head[16:24])
			width.append(w)
			height.append(h)
	
	display_path = f"{outpath}png/display/"
	if not os.path.exists(display_path) :
		print(f'Making directory:\t{display_path}')
		os.makedirs(display_path)
	
	files = [f"https://me.lsst.eu/rsolomon{file.removeprefix('/sps/lsst/users/rsolomon/web')}" for file in files]
	
	text  = f'<html>\n'
	text += f'<h2>{description}</h2>\n\n'
	
	text += f'<ul>\n'
	text += '\n'.join([f'<IMG align=center width={w} height={h} SRC={SRC}>' for w,h,SRC in zip(width,height,files)])
	text += f'\n</ul>\n\n'
	
	text += f'</html>'
	
	with open(f'{display_path}index.html', 'w') as f :
		f.write(text)

