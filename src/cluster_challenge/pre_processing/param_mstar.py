import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.modeling import models, fitting
from astropy.table import Table
import numpy as np

inpath = '/sps/lsst/users/maguena/cats/dc2/cosmoDC2_v1.1.4/mstar/'
outpath = '/sps/lsst/users/tguillem/web/clusters/mstar/plots/'

d_g = ascii.read(inpath + 'DC2_i_star.dat')
d_r = ascii.read(inpath + 'DC2_r_star.dat')
d_i = ascii.read(inpath + 'DC2_i_star.dat')
d_z = ascii.read(inpath + 'DC2_z_star.dat')
d_y = ascii.read(inpath + 'DC2_y_star.dat')

poly_g = models.Polynomial1D(degree=7)
poly_r = models.Polynomial1D(degree=7)
poly_i = models.Polynomial1D(degree=7)
poly_z = models.Polynomial1D(degree=7)
poly_y = models.Polynomial1D(degree=7)
fit = fitting.LevMarLSQFitter()

mag_star_g = fit(poly_g, d_g["col1"][1::], d_g["col2"][1::])
mag_star_r = fit(poly_r, d_r["col1"][1::], d_r["col2"][1::])
mag_star_i = fit(poly_i, d_i["col1"][1::], d_i["col2"][1::])
mag_star_z = fit(poly_z, d_z["col1"][1::], d_z["col2"][1::])
mag_star_y = fit(poly_y, d_y["col1"][1::], d_y["col2"][1::])

labels =['mag_g','mag_r','mag_i','mag_z','mag_y']
for i in range(5):
    plt.figure()
    plt.xlim(0,2.0)
    plt.ylim(14,30)
    if(i==0):
        plt.plot(d_g["col1"], d_g["col2"], 'ko', markersize=0.5)
        plt.plot(d_g["col1"], mag_star_g(d_g["col1"]))
    if(i==1):
        plt.plot(d_r["col1"], d_r["col2"], 'ko', markersize=0.5)
        plt.plot(d_r["col1"], mag_star_r(d_r["col1"]))
    if(i==2):
        plt.plot(d_i["col1"], d_i["col2"], 'ko', markersize=0.5)
        plt.plot(d_i["col1"], mag_star_i(d_i["col1"]))
    if(i==3):
        plt.plot(d_z["col1"], d_z["col2"], 'ko', markersize=0.5)
        plt.plot(d_z["col1"], mag_star_z(d_z["col1"]))
    if(i==4):
        plt.plot(d_y["col1"], d_y["col2"], 'ko', markersize=0.5)
        plt.plot(d_y["col1"], mag_star_y(d_y["col1"]))
        plt.xlabel('redshift')
    plt.ylabel(labels[i])
    plt.savefig(outpath+"m_star_"+labels[i]+".png", bbox_inches='tight')

#store parameters
coeff_g = []
coeff_g.append(mag_star_g.c0.value)
coeff_g.append(mag_star_g.c1.value)
coeff_g.append(mag_star_g.c2.value)
coeff_g.append(mag_star_g.c3.value)
coeff_g.append(mag_star_g.c4.value)
coeff_g.append(mag_star_g.c5.value)
coeff_g.append(mag_star_g.c6.value)
coeff_g.append(mag_star_g.c7.value)

coeff_r = []
coeff_r.append(mag_star_r.c0.value)
coeff_r.append(mag_star_r.c1.value)
coeff_r.append(mag_star_r.c2.value)
coeff_r.append(mag_star_r.c3.value)
coeff_r.append(mag_star_r.c4.value)
coeff_r.append(mag_star_r.c5.value)
coeff_r.append(mag_star_r.c6.value)
coeff_r.append(mag_star_r.c7.value)

coeff_i = []
coeff_i.append(mag_star_i.c0.value)
coeff_i.append(mag_star_i.c1.value)
coeff_i.append(mag_star_i.c2.value)
coeff_i.append(mag_star_i.c3.value)
coeff_i.append(mag_star_i.c4.value)
coeff_i.append(mag_star_i.c5.value)
coeff_i.append(mag_star_i.c6.value)
coeff_i.append(mag_star_i.c7.value)

coeff_z = []
coeff_z.append(mag_star_z.c0.value)
coeff_z.append(mag_star_z.c1.value)
coeff_z.append(mag_star_z.c2.value)
coeff_z.append(mag_star_z.c3.value)
coeff_z.append(mag_star_z.c4.value)
coeff_z.append(mag_star_z.c5.value)
coeff_z.append(mag_star_z.c6.value)
coeff_z.append(mag_star_z.c7.value)

coeff_y = []
coeff_y.append(mag_star_y.c0.value)
coeff_y.append(mag_star_y.c1.value)
coeff_y.append(mag_star_y.c2.value)
coeff_y.append(mag_star_y.c3.value)
coeff_y.append(mag_star_y.c4.value)
coeff_y.append(mag_star_y.c5.value)
coeff_y.append(mag_star_y.c6.value)
coeff_y.append(mag_star_y.c7.value)

t_param = Table([coeff_g,coeff_r,coeff_i,coeff_z,coeff_y], names=('mag_star_g','mag_star_r','mag_star_i','mag_star_z','mag_star_y'), meta={'name': 'mag_star'})
print(t_param)
t_param.write(outpath+'/m_star.fits', overwrite=True)
