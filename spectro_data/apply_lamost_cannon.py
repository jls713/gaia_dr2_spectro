from lamost_cannon import *
import sys
sys.path.append('/data/jls/cyanide/comparisons/')
from neural_network import *

model = tc.CannonModel.read(
    '/data/jls/GaiaDR2/spectro/lamost_cannon/lamost.cannon')
data = full_sample()#[:100000].reset_index(drop=True)
for i, l in enumerate(model.vectorizer.label_names):
    if l in data.columns:
        data[l + '_orig'] = data[l]
grid = np.load('lamost_wavelengths.npy')
output_file = '/data/jls/GaiaDR2/spectro/lamost_cannon/LAMOST_results.hdf5'

fl_dict = {}
flux = np.zeros((len(data), len(grid)))
ivar = np.zeros_like(flux)
fl_dict['flux'], fl_dict['ivar'] = flux, ivar


def load(i):
    s = load_spectrum_row(data.iloc[i])
    if s is None:
        return np.nan*np.ones(len(grid)), np.nan*np.ones(len(grid))
    return s['flux'].values, s['ivar'].values


def label_model(flux, ivar):
    if len(flux)!=len(grid) or len(ivar)!=len(grid):
        return np.nan*np.ones(7), np.nan*np.ones((7,7)), {'r_chi_sq':np.nan}
    results, cov, meta = model.test(flux, ivar, threads=44)

    ll = False
    for i in range(len(cov)):
        if cov[i] is None:
            shp = np.shape(results)[1]
            cov[i] = np.ones((shp, shp)) * np.nan
            ll = True
    if(ll):
        cov = np.stack(cov.ravel())
    return results, cov, meta


p = Pool(44)
flux, ivar = zip(*p.map(load, np.arange(len(data))))
p.close()
p.join()

snr = np.nanmedian(np.sqrt(ivar) * flux, axis=1)
data['snr']=snr
snr_sc = np.load('lamost_cannon_params.dat.npy')

result, cov, meta = label_model(flux, ivar)
for i, l in enumerate(model.vectorizer.label_names):
    data[l] = result.T[i]
    data[l + '_ERR'] = np.sqrt(cov[:, i, i])
    data[l + '_ERR_SNR'] = snr_sc[i][0] / (np.power(snr,snr_sc[i][1])+snr_sc[i][2]**2)
data['r_chi_sq'] = [m['r_chi_sq'] for m in meta]
data['snr_Cannon'] = snr
#print 'Calculate in convex hull'
#data['in_convex_hull'] = model.in_convex_hull(
#    data[list(model.vectorizer.label_names)])
print 'Loading NN'
with open('/data/jls/cyanide/comparisons/neural_network.pkl', 'r') as f:
    nn = pickle.load(f)
flds = ['C_M', 'N_M', 'ALPHA_M', 'TEFF', 'LOGG', 'M_H']
inputData = data[flds].values
inputErrData = data[[fld + '_ERR_SNR' for fld in flds]].values
results = nn.label(inputData, inputErrData, samples=100)
data['mass'], data['mass_error'] = results

data = data[['obsid','mass','mass_error','r_chi_sq','snr']+flds+[f+'_ERR_SNR' for f in flds]]

data.to_hdf(output_file, 'data')
