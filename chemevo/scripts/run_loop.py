import json
import subprocess

def outflow_radialflow():
	paramsFile = json.load(open('run_params_default.json','r'))
	outflow = [0.,0.05,0.1,0.2,0.4]
	alpha = [0.02,0.05,0.1,0.15,0.2]

	for o in outflow:
		for a in alpha:
			paramsFile['flows']['radialflow']['PezzulliAlpha']=a
			paramsFile['flows']['outflow']['Feject_in']=o
			paramsFile['flows']['outflow']['Feject_out']=o
			paramsFile['elements']=['Fe']
			with open('run_params.json','w') as outfile:
				json.dump(paramsFile,outfile)
			subprocess.call(['./run.exe','/data/jls/chem_evo/output/%iO_%iA.h5'%(int(o*100.),int(a*100.))])

def imf():
	paramsFile = json.load(open('run_params_default.json','r'))
	imf = ['Kroupa','Chabrier','Scalo','Salpeter','Tinsley']

	for i in imf:
		paramsFile['fundamentals']['IMF']=i
		o,a=0.4,0.2
		paramsFile['flows']['radialflow']['PezzulliAlpha']=a
		paramsFile['flows']['outflow']['Feject_in']=o
		with open('run_params.json','w') as outfile:
			json.dump(paramsFile,outfile)
		subprocess.call(['./run.exe','/data/jls/chem_evo/output/IMF_%s_%iO_%iA.h5'%(i,int(o*100.),int(a*100.))])

if __name__ == '__main__':
	imf()
