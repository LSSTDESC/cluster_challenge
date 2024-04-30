
import numpy as np
from scipy.special import erf



## Default parameters
r0 = 35		## richness norm.
z0 = 0.8	## redshift norm.


def log_errors(x, xerr) :
	return xerr / x / np.log(10)

def gauss(x, a, mu, sigma) :
        return a / (sigma*np.sqrt(2*np.pi)) * np.exp(-(x - mu)**2 / (2*sigma**2))

def bimodal_gauss(x, a1, mu1, sigma1, a2, mu2, sigma2) :
        return gauss(x, a1, mu1, sigma1) + gauss(x, a2, mu2, sigma2)

def power(x, a, b) :
        return a * np.power(x, b)

def mass_richness(logm0, logr, f, rnorm=r0) :
        return logm0 + f * (logr - np.log10(rnorm))

def mass_richness_zdep(logm0, logr, z, f, g, rnorm=r0, znorm=z0) :
        return logm0 + f * (logr - np.log10(rnorm)) + g * np.log10((1+z) / (1+znorm))

def polynomial(x, *B) :
	return sum([B[i] * x**i for i in range(len(B))])

def linear(x, m, b) :
	return m*x + b

def lawnchair(x, c0, c1, c2, c3) :
	## _/
	x_ = (x - c2) / c1
	return 0.5 * c0 * (x_ * erf(x_) + np.exp(-x_**2) + x_) + c3


def fit(x, y, xerr=None, yerr=None, method='least squares') :
	
	
	if xerr is None :
		xerr = np.ones_like(x)
	if yerr is None :
		yerr = np.ones_like(y)

	
