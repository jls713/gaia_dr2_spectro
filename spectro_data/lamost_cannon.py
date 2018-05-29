import urllib
from astropy.io import fits
import pandas as pd
import numpy as np
import sys
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
from apogee import load_data as load_apogee_data
sys.path.append('/data/jls/cyanide/')
from spec_utils import rebin_spectrum
from cannon_ages import *
from multiprocessing import Pool
from functools import partial


def download_spectrum(planid, lmjd, spid, fiberid):
    spectrum = 'spec-%s-%s_sp%02d-%03d.fits.gz' % (lmjd, planid, spid, fiberid)
    lamost_server = 'http://dr3.lamost.org/sas/fits/'
    fldr = '%s/' % planid
    local_location = '/media/Seagate Expansion Drive/dr3_spectra/'
    if not os.path.exists(local_location + fldr):
        os.makedirs(local_location + fldr)
    print local_location + fldr + spectrum
    if not os.path.exists(local_location + fldr + spectrum):
        urllib.urlretrieve(lamost_server + fldr + spectrum,
                           local_location + fldr + spectrum)
    return


def download_spectrum_row(row):
    ''' Pass in d.iloc[idx] '''
    return download_spectrum(row['planid'],
                             row['lmjd'],
                             row['spid'],
                             row['fiberid'])


def download_spectra_loop(lb):
    for i in range(len(lb)):
        print i
        download_spectrum_row(lb.loc[i])


def download(lb, i):
    print i
    download_spectrum_row(lb.loc[i])


def download_spectra_loop_parallel(lb, nthreads=16):
    pool = Pool(nthreads)
    pool.map(partial(download, lb), range(len(lb)))


def normalize_spectrum(spectrum, sig=50.):  # Angstrom
    c2d = np.exp(-(spectrum['wav'][:, np.newaxis] - spectrum['wav']
                   [np.newaxis, :])**2 / sig**2)
    snorm = np.dot(c2d, spectrum['flux'] * spectrum['ivar']) /\
        np.dot(c2d, spectrum['ivar'])
    spectrum['flux'] /= snorm
    spectrum['ivar'] *= snorm**2
    return spectrum


def load_spectrum(planid, lmjd, spid, fiberid):
    spectrum = 'spec-%s-%s_sp%02d-%03d.fits.gz' % (lmjd, planid, spid, fiberid)
    local_location = '/media/Seagate Expansion Drive/dr3_spectra/'
    fldr = '%s/' % planid
    print local_location + fldr + spectrum

    local_norm_location = '/data/jls/lamost/dr3_spectra/'
    if os.path.exists((local_norm_location + fldr + spectrum)[:-7] + 'hdf5'):
        return pd.read_hdf(
            (local_norm_location + fldr + spectrum)[:-7] + 'hdf5')
    else:
        return None


    hdr = fits.open(local_location + fldr + spectrum)[0].header
    data = fits.open(local_location + fldr + spectrum)[0].data.T
    data = pd.DataFrame(data.byteswap().newbyteorder(),
                        columns=['flux', 'ivar', 'wav', 'andmask', 'ormask'])
    data.loc[data['andmask'] > 0, 'ivar'] = 0.

    # Redshift
    data['wav'] = data['wav'] / (1. + hdr['Z'])

    # Rebin
    grid = np.load('lamost_wavelengths.npy')
    wav, flux, ivar = rebin_spectrum(
        data['wav'].values, data['flux'].values, grid, data['ivar'].values)
    data = pd.DataFrame(np.vstack((flux, ivar, wav)).T,
                        columns=['flux', 'ivar', 'wav'])
    n = normalize_spectrum(data)
    if not os.path.exists(local_norm_location + fldr):
        os.makedirs(local_norm_location + fldr)
    n[['flux', 'ivar']].to_hdf(
        (local_norm_location + fldr + spectrum)[:-7] + 'hdf5', 'data',
        complevel=9)
    return n


def load_spectrum_row(row):
    ''' Pass in d.iloc[idx] '''
    return load_spectrum(row['planid'],
                         row['lmjd'],
                         row['spid'],
                         row['fiberid'])

# Form dataset


def normalize_spectra_loop(lb):
    for i in range(len(lb)):
        print i
        load_spectrum_row(lb.loc[i])


def normalize(lb, i):
    print i
    load_spectrum_row(lb.loc[i])
    return


def normalize_spectra_loop_parallel(lb, nthreads=16):
    pool = Pool(nthreads)
    pool.map(partial(normalize, lb), range(len(lb)))


def generate_parent_sample():
    lamost = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_input.hdf5')
    apogee = load_apogee_data()
    fltr = apogee['ASPCAPFLAG'] == 0
    # flds = ['TEFF', 'LOGG', 'M_H', 'ALPHA_M', 'C_M', 'N_M']
    fltr &= apogee['C_M_FLAG'] == 0
    fltr &= apogee['N_M_FLAG'] == 0
    fltr &= (apogee['AK'] == apogee['AK']) & (apogee['AK'] > -1.)
    apogee = apogee[fltr].reset_index(drop=True)
    ap_c = SkyCoord(ra=apogee.RA.values * u.deg, dec=apogee.DEC.values * u.deg)
    la_c = SkyCoord(ra=lamost.ra.values * u.deg, dec=lamost.dec.values * u.deg)
    c, d2, d3 = ap_c.match_to_catalog_sky(la_c)
    apogee_cm = apogee[d2 < 2. * u.arcsec].reset_index(drop=True)
    lamost_cm = lamost.iloc[c][d2 < 2. * u.arcsec].reset_index(drop=True)
    lamost_flds = ['obsid', 'planid', 'lmjd', 'fiberid', 'spid', 'logg', 'teff']
    apogee_flds = ['APOGEE_ID', 'TEFF', 'LOGG',
                   'M_H', 'ALPHA_M', 'C_M', 'N_M', 'snr', 'AK']
    data = pd.concat((apogee_cm[apogee_flds], lamost_cm[lamost_flds]), axis=1)
    fltr = True
    for i in apogee_flds:
        fltr &= (data[i] == data[i])
    fltr &= (data['logg'] < 3.9)
    fltr &= (data['LOGG'] < 3.9)
    return data[fltr].reset_index(drop=True)


def load_spectra():
    labelled = generate_parent_sample()
    labelled.to_hdf('/data/jls/GaiaDR2/spectro/LAMOST_cannon_sample.hdf5',
                    'data')
    # return
    labelled = pd.read_hdf(
        '/data/jls/GaiaDR2/spectro/LAMOST_cannon_sample.hdf5')
    grid = np.load('lamost_wavelengths.npy')
    flux = np.zeros((len(labelled), len(grid)))
    ivar = np.zeros_like(flux)
    for i in range(len(labelled)):
        print i
        s = load_spectrum_row(labelled.iloc[i])
        flux[i, :], ivar[i, :] = s['flux'], s['ivar']
    fldr = '/data/jls/GaiaDR2/spectro/lamost_cannon/'
    labelled.to_hdf(fldr + 'training_data.hdf5', 'data')
    np.save(fldr + 'training_flux.npy', flux)
    np.save(fldr + 'training_ivar.npy', ivar)
    np.save(fldr + 'wavelengths.npy', grid)
    return


def run_cannon_model():

    fldr = '/data/jls/GaiaDR2/spectro/lamost_cannon/'
    labelled_set = pd.read_hdf(fldr + 'training_data.hdf5')
    #labelled_set = Table.from_pandas(labelled_set)
    flux = np.load(fldr + 'training_flux.npy')
    ivar = np.load(fldr + 'training_ivar.npy')
    wavelengths = np.load(fldr + 'wavelengths.npy')

    # Remove some dodgy pixels
    fltr = (flux != flux) | (ivar != ivar) | (ivar < 0.)
    flux[fltr] = 1.
    ivar[fltr] = 0.
    
    labels = ['TEFF', 'LOGG', 'M_H', 'ALPHA_M', 'C_M', 'N_M', 'AK']
    
    labelled_set, flux, ivar = pruner(labelled_set, labels, flux,
                                      ivar, wavelengths, factor=3.,
                                      threads=88)
    labelled_set.to_hdf(fldr + 'training_data_pruned.hdf5', 'data')
    np.save(fldr + 'training_flux_pruned.npy', flux)
    np.save(fldr + 'training_ivar_pruned.npy', ivar)

    test, cov, meta = \
        run_model(labelled_set, labels,
                  flux, ivar,
                  wavelengths, name='lamost',
                  output_model_file=fldr + 'lamost.cannon',
                  order=2,
                  regularization=0.,
                  photometry=None,
                  use_censors=False, threads=88,
                  generate_plots=False)
    np.save(fldr + 'test.npy', test)
    np.save(fldr + 'test_cov.npy', cov)


def full_sample():
    # But need to cut out bad fits -- upper MS
    data = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_input.hdf5')
    return data[(data.teff<6000.)&(data.logg < 3.9)].reset_index(drop=True)


if __name__ == '__main__':
    load_spectra()
    run_cannon_model()
