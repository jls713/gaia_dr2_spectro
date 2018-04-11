import aa_py
import sys
sys.path.append('../py/')
import edf_py
import json
import numpy as np
import matplotlib.pyplot as plt
try:
	import seaborn
except:
	pass

with open('../config.json') as data_file:
    data = json.load(data_file)

pot = aa_py.GalPot(str(data['potential']))
acts = aa_py.Actions_AxisymmetricFudge_InterpTables(pot,str(data['actions']),False,0.2,20.)
edf = edf_py.edf(pot,acts)
edf.readParams(str(data['edf_params']))

X = np.array([8.,2.,1.,0.1,0.1,0.1])
age = 1.
metal = -0.5

vr = np.linspace(-100.,100.,1000)
R = np.sqrt(X[0]**2+X[1]**2)
ct,st = X[0]/R, X[1]/R
vp=220.

fvr = np.array([edf(np.append(X[:3],np.array([v*ct-vp*st,v*st+vp*ct,X[-1]])),
               age,metal) for v in vr])

fvp = np.array([edf(np.append(X[:3],np.array([-(v+vp)*st,(vp+v)*ct,X[-1]])),
               age,metal) for v in vr])
fvz = np.array([edf(np.append(X[:3],np.array([-vp*st,vp*ct,v])),
               age,metal) for v in vr])

plt.plot(vr,fvr,label=r'$v_R$')
plt.plot(vr,fvp,label=r'$v_\phi$')
plt.plot(vr,fvz,label=r'$v_z$')
plt.xlabel(r'$v_r,v_\phi-220,v_z/\,\mathrm{km\,s}^{-1}$')
plt.ylabel(r'$f(v_i)$')
plt.legend()
plt.show()
