import numpy as np
import sys
sys.path.append('ezpadova')
import padova
import json

# Two age grids -- coarser for <0.5Gyr

t0_y, t1_y, dt_y = 6.6, 8.7, 0.05
t0, t1, dt = 8.725, 10.13, 0.025
Zmin, Zmax = 0.000100001, 0.06
Zsun = 0.0152
dMH = 0.01
MH = np.arange(np.log10(Zmin / Zsun), np.log10(Zmax / Zsun), dMH)
Z = np.power(10., MH) * Zsun

with open('config.json') as config:
    output_folder = json.load(config)['dir']['isochrones']+'PARSEC_Gaia_COLIBRI/'

np.savetxt(output_folder + 'Z_vals.dat', Z)
np.savetxt(output_folder + 'age_vals.dat',
           np.power(10., np.arange(t0, t1, dt)) / 1e9)

ZHgrid = np.zeros(len(Z))

phot_sys = ['sloan', 'gaiaDR2maiz', '2mass_spitzer_wise', 'ubvrijhk', 'panstarrs1']
#phot_sys = ['gaiaDR2']

for p in phot_sys:
    for n, i in enumerate(Z):
        r = padova.get_t_isochrones(
            t0, t1, dt, i, model='parsec12s_s35', ret_table=False, phot=p)
        r_y = padova.get_t_isochrones(
            t0_y, t1_y, dt_y, i, model='parsec12s_s35', ret_table=False, phot=p)
        ff = r.find('# Zini')
        gg = r[ff:].find('\n')
        ZHgrid[n] = r[ff + gg + 1:ff + gg + 1 + 6]
        psave=p
        if p is 'gaiaDR2maiz':
            psave='gaia'
        output_file = open(output_folder + 'grid/' + psave +
                           "_" + str(ZHgrid[n]) + '.dat', 'w')
        ff = r.find('# Zini')
        output_file.write(r_y+r[ff:])
        output_file.close()

    if p == 'sloan':
        np.savetxt(output_folder + 'metal_vals.dat', ZHgrid)
