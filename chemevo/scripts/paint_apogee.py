import pandas as pd
import chemevo as cc
import numpy as np

data = pd.read_csv('/data/jls/edf/apogee_sample/full_sample_Padova_we.dat')

chem_evo_data = cc.chem_evo_data('tmp.h5')

elements = ['Fe','Mg','Ca']

results = chem_evo_data.paint(data.RcP.values,data.age.values,elements)
for n,e in enumerate(elements):
	data[e]=results[n]

data.to_csv('/data/jls/edf/apogee_sample/full_sample_Padova_we.chem')
