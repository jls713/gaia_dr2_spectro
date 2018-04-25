import tabulate
import numpy as np
import glob

tablefmt = "simple"

iso_lbl = {'2MASS': r'$(J-K_s)$', 'SDSS': r'$(g-r)$',
           'Landolt': r'$(B-V)$',
           'Pan-STARRS': r'$(g_\mathrm{P}-r_\mathrm{P})$',
           'GBP_GRP': r'$(G_\mathrm{BP}-G_\mathrm{RP})$',
           'G_GRP': r'$(G-G_\mathrm{RP})$',
           'GBP_G': r'$(G_\mathrm{BP}-G)$',
           'Teff': r'$\log_{10}(T_\mathrm{eff})-4$'}


fl = np.array([np.float64(open(f).readlines()[1].split(',')) for f in glob.glob('*2015.poly')])
fl2 = np.array([np.float64(open(f).readlines()[2].split(',')) for f in glob.glob('*2015.poly')])
names = np.array([iso_lbl[f.split('_2')[0]] for f in glob.glob('*2015.poly')])

poly1 = np.vstack((names,fl.T)).T
poly2 = np.vstack((names,fl2.T)).T

print tabulate.tabulate(poly1, headers=['Colour']+['p_%i'%i for i in range(10)], tablefmt=tablefmt)
print tabulate.tabulate(poly2, headers=['Colour']+['p_%i'%i for i in range(4)], tablefmt=tablefmt)

fl = np.array([np.float64(open(f).readlines()[1].split(',')) for f in glob.glob('*2017.poly')])
fl2 = np.array([np.float64(open(f).readlines()[2].split(',')) for f in glob.glob('*2017.poly')])
names = np.array([iso_lbl[f.split('_2')[0]] for f in glob.glob('*2017.poly')])

poly1 = np.vstack((names,fl.T)).T
poly2 = np.vstack((names,fl2.T)).T

print tabulate.tabulate(poly1, headers=['Colour']+['p_%i'%i for i in range(10)], tablefmt=tablefmt)
print tabulate.tabulate(poly2, headers=['Colour']+['p_%i'%i for i in range(4)], tablefmt=tablefmt)

import pandas as pd
g = pd.read_csv('extinction_coeffs_2017.dat', skiprows=1, sep=r'\s+', header=None).values
g = g[np.argsort(g.T[0])]
g2 = pd.read_csv('extinction_coeffs_2017_leff.dat', skiprows=1, sep=r'\s+', header=None).values
g2 = g2[np.argsort(g2.T[0])]
A = np.vstack((g.T[0],g.T[7],g2.T[7]))
print tabulate.tabulate(A.T[np.argsort(g.T[7])[::-1]],
                        headers=['Band', 'A_i/E_unit_Green2018', 'A_i/E_unit_Green2018(lambda_eff)'])

g = pd.read_csv('extinction_coeffs_2015.dat', skiprows=1, sep=r'\s+', header=None).values
g = g[np.argsort(g.T[0])]
g2 = pd.read_csv('extinction_coeffs_2015_leff.dat', skiprows=1, sep=r'\s+', header=None).values
g2 = g2[np.argsort(g2.T[0])]
A = np.vstack((g.T[0],g.T[6],g2.T[6]))
print tabulate.tabulate(A.T[np.argsort(g.T[6])[::-1]],
                        headers=['Band', 'A_i/E(B-V)_Sch', 'A_i/E(B-V)_Sch(lambda_eff)'])
