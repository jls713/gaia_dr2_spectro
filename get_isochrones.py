import numpy as np
import sys
sys.path.append('ezpadova')
import padova

t0, t1, dt = 6.6, 10.13, 0.05
Zmin, Zmax = 0.000100001, 0.05
Zsun = 0.0152
dMH = 0.01
MH = np.arange(np.log10(Zmin / Zsun), np.log10(Zmax / Zsun), dMH)
Z = np.power(10., MH) * Zsun

with open('config.json') as config:
    output_folder = json.load(config)['dir']['isochrones']+'PARSEC_Gaia/'

np.savetxt(output_folder + 'Z_vals.dat', Z)
np.savetxt(output_folder + 'age_vals.dat',
           np.power(10., np.arange(t0, t1, dt)) / 1e9)

ZHgrid = np.zeros(len(Z))

phot_sys = ['sloan', 'gaia', '2mass_spitzer_wise', 'ubvrijhk', 'panstarrs1']
phot_sys = ['gaiaDR2']



for p in phot_sys:
    for n, i in enumerate(Z):
        r = padova.get_t_isochrones(
            t0, t1, dt, i, model='parsec12s', ret_table=False, phot=p)
        ff = r.find('[M/H]')
        ZHgrid[n] = r[ff + 8:ff + 14]
        psave=p
        if p is 'gaiaDR2':
            psave='gaia'
        output_file = open(output_folder + 'grid/' + psave +
                           "_" + str(ZHgrid[n]) + '.dat', 'w')
        output_file.write(r)
        output_file.close()

    if p == 'sloan':
        np.savetxt(output_folder + 'metal_vals.dat', ZHgrid)
