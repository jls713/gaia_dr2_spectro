from lamost_cannon import *
from scipy.optimize import leastsq

model = tc.CannonModel.read(
    '/data/jls/GaiaDR2/spectro/lamost_cannon/lamost.cannon')
fldr = '/data/jls/GaiaDR2/spectro/lamost_cannon/'
labelled_set = pd.read_hdf(fldr + 'training_data_pruned.hdf5')
labelled_set = Table.from_pandas(labelled_set)
test = np.load(fldr + 'test.npy')
cov = np.load(fldr + 'test_cov.npy')
flux = np.load(fldr + 'training_flux_pruned.npy')
ivar = np.load(fldr + 'training_ivar_pruned.npy')

snr = np.nanmedian(flux * np.sqrt(ivar), axis=1)

from plotting_general import running_median
gg = np.zeros((7,3))
for i, n in enumerate(model.vectorizer.label_names):
    #bin_edge = np.array([0., 20., 40., 60., 80., 100., 120.])
    #bins = np.digitize(snr, bin_edge)
    #bin_centres = .5 * (bin_edge[1:] + bin_edge[:-1])
    #std = [np.nanstd((test.T[i] - labelled_set[n])[bins == j])
    #     for j in range(1, len(bin_edge))]
    #stdI = [np.nanmedian(cov[:,i,i][bins == j])
    #     for j in range(1, len(bin_edge))]
    #print np.nanstd((test.T[i] - labelled_set[n]))

    r = running_median(snr, test.T[i]-labelled_set[n])
    std = .5*(r[3]-r[2])
    bin_centres = r[0]
    def fn(p):
        return std - p[0] / (np.power(bin_centres,p[1])+p[2]**2)
    gg[i]=leastsq(fn, [1.,.5,0.2])[0]
np.save('lamost_cannon_params.dat',gg)
