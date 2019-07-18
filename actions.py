import numpy as np
import sys
sys.path.append('edf_sampling/py')
import edf_sampling as edf_sampling
from multiprocessing import Pool
from functools import partial


def init_():
    edf_sampling.setup(False, False)


in_data = {}
npool = 200
Nsamples = 50


def rejection(samples, mu_old, mu_new, sig_old, sig_new):
    xmin = mu_old - 3. * sig_old
    if(mu_new - 3. * sig_new < xmin):
        xmin = mu_new - 3. * sig_new
    xmax = mu_old + 3. * sig_old
    if(mu_new + 3. * sig_new > xmax):
        xmax = mu_new + 3. * sig_new
    xx = np.linspace(xmin, xmax, 1000)
    mm = np.exp(-(xx - mu_new)**2 / 2. / sig_new**2) / \
        np.exp(-(xx - mu_old)**2 / 2. / sig_old**2)
    M = np.max(mm)
    mm = np.exp(-(samples[:, 0] - mu_new)**2 / 2. / sig_new**2) / \
        np.exp(-(samples[:, 0] - mu_old)**2 / 2. / sig_old**2)
    nn = np.random.uniform(size=len(samples))
    return samples[nn < mm / M]


def rejection_dm(samples, mu_old, mu_new, sig_old, sig_new):
    xmin = mu_old - 3. * sig_old
    if(1. / np.power(10., 0.2 * (mu_new + 3. * sig_new) - 2.) < xmin):
        xmin = 1. / np.power(10., 0.2 * (mu_new + 3. * sig_new) - 2.)
    xmax = mu_old + 3. * sig_old
    if(1. / np.power(10., 0.2 * (mu_new - 3. * sig_new) - 2.) > xmax):
        xmax = 1. / np.power(10., 0.2 * (mu_new - 3. * sig_new) - 2.)
    xx = np.linspace(xmin, xmax, 1000)
    mm = np.exp(-(5. * np.log10(100. / xx) - mu_new)**2 / 2. /
                sig_new**2) / np.exp(-(xx - mu_old)**2 / 2. / sig_old**2)
    M = np.max(mm)
    mm = np.exp(-(5. * np.log10(100. / samples[:, 0]) - mu_new)**2 / 2. /
                sig_new**2) / np.exp(-(samples[:, 0] - mu_old)**2 / 2. / sig_old**2)
    nn = np.random.uniform(size=len(samples))
    return samples[nn < mm / M]


def actions_samples(i):
    inputs, outputs = in_data['inputs'], in_data['outputs']
    if outputs['par'][i] != outputs['par'][i] or outputs['par_err'][i] != outputs['par_err'][i] or inputs['parallax'][i] != inputs['parallax'][i] or inputs['pmra_error'][i] != inputs['pmra_error'][i] or inputs['hrv'][i] != inputs['hrv'][i] or inputs['e_hrv'][i] != inputs['e_hrv'][i] or inputs['e_hrv'][i] <= 0.:
        return np.zeros(28) * np.nan
    covar = np.zeros((3, 3))
    covar[0][0] = inputs['parallax_error'][i]**2
    covar[1][1] = inputs['pmra_error'][i]**2
    covar[2][2] = inputs['pmdec_error'][i]**2
    covar[0][1] = covar[1][0] = inputs['parallax_pmra_corr'][i] * \
        inputs['pmra_error'][i] * inputs['parallax_error'][i]
    covar[0][2] = covar[2][0] = inputs['parallax_pmdec_corr'][i] * \
        inputs['pmdec_error'][i] * inputs['parallax_error'][i]
    covar[1][2] = covar[2][1] = inputs['pmra_pmdec_corr'][i] * \
        inputs['pmra_error'][i] * inputs['pmdec_error'][i]
    fr = []
    N = 1000
    while len(fr) < Nsamples and N < 200000:
        sampl = np.random.multivariate_normal(
            np.array([inputs['parallax'][i], inputs['pmra']
                      [i], inputs['pmdec'][i]]),
            covar, size=N)
        fr = rejection_dm(sampl, inputs['parallax'][i], outputs['dm'][i],
                          inputs['parallax_error'][i], outputs['dm_err'][i])
        N *= 2
    fr[:, 0] = 5. * np.log10(100. / fr[:, 0])
    if len(fr) > Nsamples:
        fr = fr[:Nsamples]
    else:
        fr = np.zeros((Nsamples, 3))
        fr[:, 0] = np.random.normal(
            outputs['dm'][i], outputs['dm_err'][i], Nsamples)
        fr[:, 1:] = sampl[:Nsamples, 1:]
    samples = np.zeros((Nsamples, 4))
    samples[:, 0] = np.power(10., 0.2 * fr[:, 0] - 2.)
    samples[:, 1:3] = fr[:, 1:]
    samples[:, 3] = np.random.normal(
        inputs['hrv'][i], inputs['e_hrv'][i], Nsamples)
    print i, np.nanmedian(samples, axis=0)
    Xs = np.array([
                  edf_sampling.process_data(
                      np.array([np.deg2rad(inputs['ra'][i]),
                                np.deg2rad(inputs['dec'][i]),
                                samples[ns][0],
                                samples[ns][3],
                                samples[ns][1],
                                samples[ns][2]]))
                  for ns in range(Nsamples)
                  ])
    if np.count_nonzero(~np.isinf(Xs[:, -1])) == 0:
        XX = np.zeros(28) * np.nan
        XX[:10]=np.nanmean(Xs,axis=0)[2:12]
        XX[14:24]=np.nanstd(Xs,axis=0)[2:12]
    else:
        XX = np.concatenate((
                             np.nanmean(Xs, axis=0)[2:12],
                             np.nanmean(Xs[~np.isinf(Xs[:, -1])], axis=0)[12:],
                             np.nanstd(Xs, axis=0)[2:12],
                             np.nanstd(Xs[~np.isinf(Xs[:, -1])], axis=0)[12:]
                           ))
    print i, XX
    return XX


additional_output = ['s',
                     'vlos', 'mu_l', 'mu_b',
                     'R', 'phi', 'z',
                     'vR', 'vphi', 'vz',
                     'JR', 'Lz', 'Jz',
                     'Rc']
additional_output_errors = [a + '_err' for a in additional_output]


def process_actions(input_data, output_data, output_path=None):

    in_data['inputs'] = input_data
    in_data['outputs'] = output_data.to_pandas()
    num_tasks = len(input_data)
    nindex = range(num_tasks)
    nsplit = 2
    for Nindx in np.array_split(nindex, nsplit):
        if(npool > 1):
            p = Pool(npool)
            r = p.map(partial(actions_samples),
                      Nindx,
                      )
            p.close()
            p.join()
        else:
            r = map(lambda i: actions_samples(
                i), Nindx)
        # for nindx, X in zip(Nindx, r):
        output_cols = additional_output + additional_output_errors
        for K, n in enumerate(output_cols):
            output_data[n][Nindx] = np.array(r)[:,K]
        if output_path:
            output_data.write(output_path, path='data', format='hdf5',
                              compression=True,
                              serialize_meta=True,
                              overwrite=True)
    print 'All actions computed!'
    return output_data
