
import numpy as np
import healpy as hp
import yaml, os, sys, shutil
import clevar
from astropy.table import Table

src = __file__.replace('mt_utils.py', '')

class mtCatalogs :
    """
    Catalog matching object for matching catalogs, members, and footprints based on a single configuration file.
    """

    def __init__(self, config) :
        if isinstance(config, dict) :
            self.cfg = config
        elif isinstance(config, str) :
            with open(config) as fstream :
                cfg = yaml.safe_load(fstream)
            self.cfg = cfg

    def get_clevar_version(self) :
        return(clevar.version.__version__)

    def get_catalogs(self) :
        self.cats = {}
        for cat in ['cat1', 'cat2'] :
            self.cats[cat] = _load_file(
                    self.cfg['paths']['matching']['in'][cat]['cluster'],
                    method = 'clevar',
                    name = self.cfg['cats'][cat],
                    tags = self.cfg['cat_keys'][cat]['cluster'],
                    full = False
                    )

    def get_mb_catalogs(self) :
        if hasattr(self, 'cats') :
            for cat in ['cat1', 'cat2'] :
                mb_tags = self.cfg['cat_keys'][cat]['member']
                mb_cat = Table.read(self.cfg['paths']['matching']['in'][cat]['member'])
                self.cats[cat].add_members(
                        id = mb_cat[mb_tags['id']].astype(str),
                        id_cluster = mb_cat[mb_tags['id_cluster']].astype(type(self.cats[cat]['id'][0])),
                        ra = mb_cat[mb_tags['ra']],
                        dec = mb_cat[mb_tags['dec']],
                        z = mb_cat[mb_tags['z']],
                        mag_u = mb_cat[mb_tags['mag_u']],
                        mag_g = mb_cat[mb_tags['mag_g']],
                        mag_r = mb_cat[mb_tags['mag_r']],
                        mag_i = mb_cat[mb_tags['mag_i']],
                        mag_z = mb_cat[mb_tags['mag_z']],
                        mag_y = mb_cat[mb_tags['mag_y']]
                        )
        else :
            raise ValueError(f"Must first read in cluster catalog with 'self.get_catalogs()' function.")


    def apply_mt_footprints_mask(self) :
        
        ftps = {}
        for cat in ['cat1', 'cat2'] :
            ftps[cat] = _load_file(os.path.join(src, f"../footprints/{self.cfg['footprints'][cat]['pixFile']}"))

        hpix = {'cat1':
                hp.ang2pix(
                    self.cfg['footprints']['cat2']['NSIDE'],
                    self.cats['cat1'][self.cfg['cat_keys']['cat1']['cluster']['ra']],
                    self.cats['cat1'][self.cfg['cat_keys']['cat1']['cluster']['dec']],
                    nest=False,
                    lonlat=True),
                'cat2':
                hp.ang2pix(
                    self.cfg['footprints']['cat1']['NSIDE'],
                    self.cats['cat2'][self.cfg['cat_keys']['cat2']['cluster']['ra']],
                    self.cats['cat2'][self.cfg['cat_keys']['cat2']['cluster']['dec']],
                    nest=False,
                    lonlat=True)
                }

        self.cats['cat1'] = self.cats['cat1'][np.isin(hpix['cat1'], ftps['cat2'])]
        self.cats['cat2'] = self.cats['cat2'][np.isin(hpix['cat2'], ftps['cat1'])]


    def match_catalogs(self) :
        if self.cfg['matching']['method'] == 'proximity' :
            mt = clevar.match.ProximityMatch()
        elif self.cfg['matching']['method'] == 'member' :
            mt = clevar.match.MembershipMatch()

        cosmo = clevar.cosmology.AstroPyCosmology()
        cosmo._init_from_params(
                H0 = self.cfg['cosmology']['H0'],
                Omega_b0 = self.cfg['cosmology']['Ob0'],
                Omega_dm0 = self.cfg['cosmology']['Om0'],
                Omega_k0 = self.cfg['cosmology']['Ok0']
                )

        mt.match_from_config(self.cats['cat1'], self.cats['cat2'], self.cfg['matching']['match_config'], cosmo=cosmo)


    def write_mtCatalogs(self) :
        outpath = _make_output_pathNames_from_config(self.cfg)
        
        if os.path.exists(outpath) :
            shutil.rmtree(outpath)
        os.makedirs(outpath)

        self.cats['cat1'].write(f"{outpath}{self.cats['cat1'].name}.fits", overwrite=True)
        self.cats['cat2'].write(f"{outpath}{self.cats['cat2'].name}.fits", overwrite=True)

    def write_mtCatalogs_mbs(self) :
        outpath = _make_output_pathNames_from_config(self.cfg)

        self.cats['cat1'].members.write(os.path.join(outpath, self.cats['cat1'].name + "_mbs.fits"), overwrite=True)
        self.cats['cat2'].members.write(os.path.join(outpath, self.cats['cat2'].name + "_mbs.fits"), overwrite=True)



def _make_output_pathNames_from_config(cfg) :
    outpath_base = cfg['paths']['matching']['out']

    if cfg['matching']['method'] == 'proximity' :
        delta_z = float(cfg['matching']['match_config']['catalog1']['delta_z'])
        mt_radius = float(cfg['matching']['match_config']['catalog1']['match_radius'])
        mt_pref = cfg['matching']['match_config']['preference']

        outpath = os.path.join(
                outpath_base,
                f'proximity_matching/deltaz_{delta_z}_matchradius_{mt_radius}mpc_pref_{mt_pref}/'
                )
    elif cfg['matching']['cl_method'] == 'member' :
        delta_z = float(cfg['matching']['match_config']['match_members_kwargs']['delta_z'])
        min_share_fract1 = float(cfg['matching']['match_config']['minimum_share_fraction1'])
        min_share_fract2 = float(cfg['matching']['match_fonfig']['minimum_share_fraction2'])
        mt_pref = cfg['matching']['match_config']['preference']
        mt_method_mb = cfg['matching']['match_config']['match_members_kwargs']['method']

        outpath = os.path.join(
                outpath_base,
                f'member_matching/fshare1_{min_share_fract1}_fshare2_{min_share_fract2}_pref_{mt_pref}__mbmt_{mt_method_mb}/'
                )

    return outpath


def _load_file(fileName, method='suffix', **kwargs) :
    if os.path.exists(fileName) :
        if method == 'suffix' :
            fileType = fileName.split('.')[-1]

            if (fileType == 'npy') :
                return np.load(fileName)
        elif method == 'clevar' :
            return clevar.catalog.ClCatalog.read(fileName, **kwargs)
    else :
        sys.exit(f"The file '{fileName}' does not exist.")
