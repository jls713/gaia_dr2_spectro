import pandas as pd
import chemevo as cc
import numpy as np

data = pd.read_csv('/data/jls/edf/apogee_sample/full_sample_Padova_we.dat')

chem_evo_data = cc.chem_evo_data('/data/jls/chem_evo/results/singlezone_multizone_dump.hdf5')

elements = ['Fe','Si']

results = chem_evo_data.paint(data.RcP.values,data.age.values,elements)
for n,e in enumerate(elements):
	data[e]=results[n]

data.to_csv('/data/jls/edf/apogee_sample/full_sample_Padova_we_GASDUMP.chem')
