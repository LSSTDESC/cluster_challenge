import numpy as np
import glob, os, shutil, pickle, struct

import matplotlib.pyplot as plt




def make_path(cfg, base, addon=None) :
	mt_method = cfg['matching']['method']
	mt_pref = cfg['matching']['pref']

	if mt_method == 'member' :
		mt_params = [float(cfg['matching']['minimum_share_fraction']),]
		path = os.path.join(base, 'member_matching', f"fshare_{mt_params[0]}_pref_{mt_pref}")
	elif mt_method == 'proximity' :
		mt_params = [float(cfg['matching']['delta_z']), float(cfg['matching']['match_radius']),]
		path = os.path.join(base, 'proximity_matching',
				f"deltaz_{mt_params[0]}_matchradius_{mt_params[1]}mpc_pref_{mt_pref}")

	if addon is None :
		return path
	else :
		path = os.path.join(path, addon)
		return path




def save_figure(fig, outpath, saveas, png=True, pdf=False, pkl=False) :
	ftypes = ['png', 'pdf', 'pkl']
	outpaths = {}
	for (i,ftypeBool) in enumerate([png, pdf, pkl]) :
		outpaths[ftypes[i]] = os.path.join(outpath, ftypes[i])
		if ftypeBool & (not os.path.exists(outpaths[ftypes[i]])) :
			print(f"Making directory:\t{outpaths[ftypes[i]]}")
			os.makedirs(outpaths[ftypes[i]])

	if png :
		fig.savefig(os.path.join(outpaths['png'], f"{saveas}.png"), bbox_inches='tight')
	if pdf :
		fig.savefig(os.path.join(outpaths['pdf'], f"{saveas}.pdf"), bbox_inches='tight')
	if pkl :
		pickle.dump(fig, open(os.path.join(outpaths['pkl'], f"{saveas}.pickle"), 'wb'))
	plt.close(fig)
	
	if png :
		update_index_file_for_html_display(outpath)



def update_index_file_for_html_display(outpath, description='') :
	files = list(filter(os.path.isfile, glob.glob(f'{outpath}/png/*')))
	files.sort(key=lambda x: os.path.getmtime(x))
	
	width = []
	height = []
	for f in files :
		with open(f, 'rb') as fhandle :
			head = fhandle.read(24)
			w, h = struct.unpack('>ii', head[16:24])
			width.append(w)
			height.append(h)
	
	display_path = f"{outpath}/png/display/"
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

