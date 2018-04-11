## Script for taking the data file from http://iopscience.iop.org/article/10.1086/392523/fulltext/ and processing into a simpler form to read into the chemevo code.

import numpy as np
import pandas as pd

datafile = '/data/jls/chem_evo/yields/type_II/chieffi_limongi.out'
outfile = '/data/jls/chem_evo/yields/type_II/chieffi_limongi.dat'

Z = {0.:0,1e-6:1,1e-4:2,1e-3:3,6e-3:4,2e-2:5}
M = {13.:0,15.:1,20.:2,25.:3,30.:4,35.:5}

data_frame = pd.DataFrame(0.,index=np.arange(len(M)*len(Z)),
                          columns=['M','Z','M_eject','M_remnant'])
data_frame['M'] = np.ravel([np.sort(M.keys())]*len(Z))
data_frame['Z'] = np.ravel([np.sort(Z.keys()*len(M))])

with open(datafile) as f:
	for i in range(17):
		f.readline()
	for l in f:
		model_descr = l[:33].rstrip()
		if(model_descr[0]=='R'):
			continue
		z = float(model_descr.split('=')[-1])
		model_param = l[34:44].rstrip()
		if model_param[0]=='M':
			if(model_param=='M_ejected_'):
				key='M_eject'
			if(model_param=='M_remnant_'):
				key='M_remnant'
		else:
			key = model_param.split('^')[-1]
		if key not in data_frame.columns:
			data_frame[key]=np.zeros(len(data_frame))
		for m in range(6):
			data_frame.loc[Z[z]*len(M)+m,key]+=np.float64(l[45+10*m:54+10*m].rstrip())

cols = list(data_frame)
cols.insert(cols.index('Co'),cols.pop(cols.index('Ni')))
data_frame = data_frame.ix[:,cols]
data_frame.to_csv(outfile,sep=' ',index=False)
