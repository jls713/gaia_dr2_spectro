## Script for taking the data file from http://iopscience.iop.org/article/10.1086/508914/fulltext/63689.tables.html and processing into a simpler form to read into the chemevo code.

import numpy as np
import pandas as pd

data = pd.read_csv('/data/jls/chem_evo/yields/type_II/kobayashi2006_new.dat',sep=r'\s+',names=['Z','type','13.','15.','18.','20.','25.','30.','40.'],skiprows=24)
by_z = data.groupby('Z')
f = open('/data/jls/chem_evo/yields/type_II/kobayashi2006_proc.dat', 'w').close()
hdr=True
for n,g in by_z:
    g=g.transpose()
    g.columns = g.values[1]
    g['Z']=np.ones(len(g))*g.values[0][0]
    Z = g['Z']
    g.drop(labels=['Z'], axis=1,inplace = True)
    g.insert(0, 'Z', Z)
    g=g.drop(g.index[[0,1]])
    for t in g.columns:
        if '^' in t:
            el = t.split('^')[-1]
            if el in g.keys():
                g[el]+=g[t]
            else:
                g[el]=g[t]
    for t in g.columns:
        if '^' in t:
            del g[t]
    g.insert(3, 'H', g['p']+g['d'])
    del g['p']
    del g['d']
    with open('/data/jls/chem_evo/yields/type_II/kobayashi2006_proc.dat', 'a') as f:
        g.to_csv(f,sep=' ',header=hdr)
    hdr=False
