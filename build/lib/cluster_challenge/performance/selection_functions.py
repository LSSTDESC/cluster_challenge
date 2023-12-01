# Selection functions for cosmology with cluster count on cosmoDC2/DC2 data
# Author: T. Guillemin
# Date: July 2022

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

###fit functions
def gaussian(x, mu, sigma):
    return 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))

def mu_lnlambda(redshift, logM, mu0, G_z_mu, G_logM_mu):
    #return mu0 + G_z_mu*(1+redshift)/(1+z_0) + G_logM_mu*(logM/logM_0)**0.5
    return mu0 + G_z_mu*np.log10((1+redshift)/(1+z_0)) + G_logM_mu*(logM-logM_0)

def sigma_lnlambda(redshift, logM, sigma0, F_z_sigma, F_logM_sigma):
    #return sigma0 + F_z_sigma*(1+redshift)/(1+z_0) + F_logM_sigma*(logM)/logM_0
    return sigma0 + F_z_sigma*np.log((1+redshift)/(1+z_0)) + F_logM_sigma*(logM-logM_0)

###store parameters
###Cosmo DC2
#case 1b: true richness, true redshift
#from CLCosmo_Sim package:
#def mu_loglambda_logM_f(self, redshift, logm): ... return loglambda0 + A_z_mu * np.log10((1+redshift)/(1 + z0)) + A_logm_mu * (logm-np.log10(m0))
#def sigma_loglambda_logm_f(self, redshift, logm): ... return sigma_lambda0 + A_z_sigma * np.log10((1+redshift)/(1 + z0)) + A_logm_sigma * (logm-np.log10(m0))
#Skysim_image / mass = m200c / halos with m200c>10**13
[z_0,logM_0]=[0.78,13.27]
[loglambda0,A_z_mu,A_logm_mu,sigma_lambda0,A_z_sigma,A_logm_sigma]=[2.830,1.439,1.935,0.424,-0.136,-0.127]
[err_loglambda0,err_A_z_mu,err_A_logm_mu,err_sigma_lambda0,err_A_z_sigma,err_A_logm_sigma]=[0.0009,0.013,0.003,0.0006,0.004,0.002]

###validation plots
do_plots = False
if(do_plots==False):
    sys.exit()    
outpath = "/sps/lsst/users/tguillem/web/clusters/cluster_challenge/selection_function/validation/case_1b/"
if os.path.exists(outpath):
         shutil.rmtree(outpath)
os.makedirs(outpath) 
#richness plots in (m,z) bins
zbins = [0,0.5,0.75,1.0,1.2]
ybins = np.linspace(13, 14.6, 9)
n_z = len(zbins)-1
n_y = len(ybins)-1
a_z = np.zeros((n_z,n_y))
a_mass = np.zeros((n_z,n_y))
for i in range(0,n_z):
    cut1 = zbins[i]
    cut2 = zbins[i+1]
    for j in range(0,n_y):
        cut3 = ybins[j]
        cut4 = ybins[j+1]
        a_z[i][j] = (cut1+cut2)/2
        a_mass[i][j] = (cut3+cut4)/2
        plt.figure()
        plt.xlabel("ln(r)");
        plt.ylabel("P(ln(r)|m,z)")
        f_cut1=round(cut1,1)
        f_cut2=round(cut2,1)
        f_cut3=round(cut3,1)
        f_cut4=round(cut4,1)
        plt.title('cosmoDC2: log(m200c) '+str(f_cut3)+'-'+str(f_cut4) + ' / z '+str(f_cut1)+'-'+str(f_cut2))
        #plot parametrization
        x = np.linspace(0.1, 6, 2000)
        mu = mu_lnlambda(a_z[i][j],a_mass[i][j],loglambda0,A_z_mu,A_logm_mu)
        sigma = sigma_lnlambda(a_z[i][j],a_mass[i][j],sigma_lambda0,A_z_sigma,A_logm_sigma)
        gauss_2 = gaussian(x, mu, sigma)
        plt.plot(x, gauss_2, color='red', linewidth=2.0,label="Param")
        plt.legend()
        plt.savefig(outpath+'richness_redshift_bin_'+str(i)+'_mass_bin_'+str(j)+'.png')
        plt.close()
