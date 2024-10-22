
## Used to set consistent fontsize and style across all plotting scripts.

import matplotlib as mpl



def figsize(n=1, m=1) :
        return (n*5, m*5)

def set_style(fontsize=12) :
	mpl.rcParams['font.size'] = fontsize
