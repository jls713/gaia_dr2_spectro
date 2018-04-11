import sys
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from multiprocessing import Pool
from functools import partial
from reddening import ReddeningMaps
from astropy.table import Table
import json
redd = ReddeningMaps(with_extinction_maps=True)
sys.path.append('isochrone')
sys.path.append('edf_sampling/py')
import isodist_js as isodist
import edf_sampling as edf_sampling


def init(types="All", wemap=False, mean_feh_err=1e-5):
    isodist.init_isochrone(types, 1, mean_feh_err)
    if(wemap):
        isodist.load_emap()
    isodist.load_prior("2018")
    edf_sampling.setup(False, False)


def process_single_distance(input_data, w, prior, with_parallax, with_mass, i):
    ''' Computes a single distance to the star in entry i in the table
        input_data. w is a string which gives the isochrone set to use,
        if prior we use the prior, if with_parallax we use the parallax,
        if with_mass we use the mass'''

    #  First photometric inputs
    mag_str = input_data['mag_use'][i]
    mags = np.array([np.float64(input_data[k][i]) for k in mag_str])
    errmags = np.array([np.float64(input_data['e' + k][i]) for k in mag_str])

    c_icrs = SkyCoord(ra=input_data.ra[i] * u.degree,
                      dec=input_data.dec[i] * u.degree, frame='icrs')
    ll = c_icrs.galactic.l.radian
    b = c_icrs.galactic.b.radian
    ldeg = c_icrs.galactic.l.degree
    bdeg = c_icrs.galactic.b.degree

    #  Now spectroscopic inputs
    teff, logg, Z = \
        np.log10(input_data.teff[i]), input_data.logg[i], input_data.fe_h[i]

    data = np.array([np.float64(Z), np.float64(teff), np.float64(logg),
                     np.float64(ll), np.float64(b)])

    median_log_av, std_log_av, dist_ebv = redd.get_log_av_along_los(ldeg, bdeg)

    ERRZ = input_data.e_fe_h[i]
    ERRteff = input_data.e_teff[i] / input_data.teff[i] / np.log(10.)
    ERRlogg = input_data.e_logg[i]

    CZT = input_data.rho_TZ[i] * ERRteff * ERRZ
    CGZ = input_data.rho_gZ[i] * ERRlogg * ERRZ
    CTG = input_data.rho_Tg[i] * ERRteff * ERRlogg

    data_errs = np.array([np.float64(ERRZ)**2,
                          np.float64(CZT), np.float64(CGZ),
                          np.float64(ERRteff)**2,
                          np.float64(CTG), np.float64(ERRlogg)**2])

    if(with_parallax):
        parallax = np.float64(input_data.parallax[i])
        parallax_error = np.float64(input_data.parallax_error[i])
    else:
        parallax = 0.
        parallax_error = -1.

    if(with_mass):
        mass = np.float64(input_data.mass[i])
        mass_error = np.float64(input_data.mass_error[i])
    else:
        mass = 0.
        mass_error = -1.

    ir_problem = np.sum([(errmags[k] < 0. or mags[k] < -10. or
                          mags[k] != mags[k] or errmags[k] != errmags[k])
                         for k in range(len(mags))])
    spec_problem = (logg < -5000. or logg != logg or
                    teff < -5000. or teff != teff or
                    Z <= -5000. or Z != Z or
                    ERRteff != ERRteff or
                    ERRlogg != ERRlogg or
                    ERRZ != ERRZ)
    astrom_problem = (parallax_error != parallax_error or
                      parallax != parallax or
                      (parallax_error < 0. and parallax > 0.))
    mass_problem = (mass_error != mass_error or
                    mass != mass or
                    (mass_error < 0. and mass > 0.))

    if(not spec_problem and not ir_problem and not astrom_problem):
        print i, mags, data, errmags, data_errs, mag_str
        func = isodist.prob_distance_extinctprior
        distances = func(mags, data,
                         errmags, data_errs, prior,
                         np.float64(median_log_av), np.float64(std_log_av),
                         np.log(np.float64(dist_ebv.value)),
                         w, mag_str, False, parallax, parallax_error,
                         mass, mass_error)
        # distances = np.append(distances[:-5], distances[-1])
        # Convert log() to log10()
        distances[7:9] /= np.log(10.)
        distances[13:15] /= np.log(10.)
        print i, distances
        input_data.loc[i, 's'] = 10**(0.2 * distances[1] - 2.)
        if(distances[2] < 0.):
            X = np.ones(17) * np.nan
            distances[:2] = np.nan
            distances[3:] = np.nan
            X[-1] = 1
        else:
            X = edf_sampling.process_data(
                np.array([np.deg2rad(input_data['ra'][i]),
                          np.deg2rad(input_data['dec'][i]),
                          input_data['s'][i],
                          input_data['hrv'][i],
                          input_data['pmra'][i],
                          input_data['pmdec'][i]]))
            X = np.append(X, [0])
            # X = np.delete(X, 3)  # remove radial velocity from output
        return np.append(distances[1:], X)
    else:
        XX = np.ones(36) * np.nan
        if spec_problem:
            XX[-1] = 2
        if ir_problem:
            XX[-1] = 3
        if astrom_problem:
            XX[-1] = 4
        if mass_problem:
            XX[-1] = 5
        return XX


additional_output = ['l', 'b', 's',
                     'vlos', 'mu_l', 'mu_b',
                     'R', 'phi', 'z',
                     'vR', 'vphi', 'vz',
                     'JR', 'Lz', 'Jz',
                     'Rc', 'flag']
pos_unit = u.kpc
vel_unit = u.km / u.s
action_unit = u.kpc * u.km / u.s
log10_mag_unit = u.dimensionless_unscaled

additional_output_units = [u.rad, u.rad, u.kpc,
                           vel_unit, u.mas / u.yr, u.mas / u.yr,
                           pos_unit, pos_unit, pos_unit,
                           vel_unit, vel_unit, vel_unit,
                           action_unit, action_unit, action_unit,
                           pos_unit,
                           u.dimensionless_unscaled
                           ]

distance_output = ['dm', 'dm_err',
                   'dist', 'dist_err',
                   'par', 'par_err',
                   'log10_age', 'log10_age_err',
                   'mass', 'mass_err',
                   'Z', 'Z_err',
                   'log10_av', 'log10_av_err',
                   'log10_teff', 'log10_teff_err',
                   'logg', 'logg_err',
                   'covar_dm_logage']

distance_units = [u.mag, u.mag,
                  pos_unit, pos_unit,
                  u.mas, u.mas,
                  u.dex(u.Gyr), u.dex(u.Gyr),
                  u.Msun, u.Msun,
                  u.dex, u.dex,
                  log10_mag_unit, log10_mag_unit,
                  u.dex(u.K), u.dex(u.K),
                  u.dex(u.cm / u.s**2), u.dex(u.cm / u.s**2),
                  u.dimensionless_unscaled]


def new_add_columns(data):
    for i, unit in zip(additional_output, additional_output_units):
        data[i] = np.nan * unit
    data['flag'] = data['flag'].astype(int)
    return data


def new_distance_columns(data):
    for i, unit in zip(distance_output, distance_units):
        data[i] = np.nan * unit
    return data


def column_descr(name, id_col, which, prior, with_parallax, with_mass):

    with open('config.json') as data_file:
        data = json.load(data_file)
    R0 = data["solar_motion"][0]
    z0 = data["solar_motion"][1]
    U0 = data["solar_motion"][2]
    V0 = data["solar_motion"][3]
    W0 = data["solar_motion"][4]
    pot = data["potential"].split('/')[-1]
    which_dict = {'Padova': 'PARSEC version 1.2S eta=0.2 '
                            'http://stev.oapd.inaf.it/cgi-bin/cmd',
                  'Dartmouth': 'Dartmouth 2012 '
                               'http://stellar.dartmouth.edu/models/grid.html',
                  'BaSTI': 'BaSTI canonical eta=0.4 '
                           'http://basti.oa-teramo.inaf.it/'}
    if isinstance(which, list):
        which_string = '\n'.join(['%s' % which_dict[w] for w in which])
    else:
        which_string = '%s' % which_dict[which]
    parallax_string = ''
    if with_parallax:
        parallax_string = ' combined with Gaia DR2'
    if with_mass:
        mass_string = ' using spectroscopic mass estimates'
    prior_string = 'No prior'
    if prior:
        prior_string = '2018 prior described in Sanders et al. (2018)'
    descr = ['Distances, ages, masses and extinctions for the spectroscopic'
             ' survey %s' % name + parallax_string + mass_string,
             'Isochrone set(s) used = ' + which_string,
             'Adopted prior: %s' % prior_string,
             '%s is the unique identifier for the %s survey' % (id_col, name),
             '(l,b) are Galactic coordinates (in rad)',
             's is distance (kpc) used for the Galactocentric velocity and '
             'action calculations = 10**(0.2*dm-2)',
             'vlos is radial velocity (km/s)',
             '(mu_l, mu_b) = Proper motion in Galactic coordinates (mas/yr)'
             ' (mu_l contains cos(b) factor)',
             '(R,phi,z) are Galactocentric coordinates (kpc) -- convention is'
             ' right-handed coordinate system with',
             ' positive x away from GC and positive z towards NGP',
             '(vR,vphi,vz) are Galactocentric polar velocities (km/s)-- note '
             'vphi is in the negative phi direction ',
             '(so Sun has positive vphi)',
             'We assume the Sun is at (R,z) = (%0.1f,%0.2f) kpc' % (R0, z0) +
             ' with peculiar velocity (U,V,W) = '
             '(%0.2f,%0.2f,%0.2f) km/s' % (U0, V0, W0),
             'We use the potential %s from McMillan (2017).' % pot,
             'The actions (JR,Lz,Jz) are computed using the Staeckel fudge'
             ' from Binney (2012) (in units kpc km/s).',
             'Rc is the radius of a circular orbit with angular momentum Lz '
             '(in kpc).',
             'dm and dm_err give the distance modulus with associated error',
             'dist and dist_err give the distance (kpc) with associated error',
             'par and par_err give the parallax (mas) with associated error',
             'mass and mass_err give the mass (Msun) with associated error',
             'Z and Z_err give the metallicity estimates',
             'log10_av and log10_av_err give log10 V-band extinction'
             ' This is the one column without units as ',
             ' astropy cannot handle'
             ' log10(mag)',
             'log10_teff and log10_teff_err give the log10 effective temp (K)',
             'logg and logg_err give the log surface gravity (cm/s^2)',
             'covar_dm_logage gives the covariance between DM and log10 age',
             'flag is an integer where 0 denotes all columns contain values',
             'if flag!=0 then many other columns will be nan',
             'if flag=1 the pipeline failed (could not find overlap with any'
             ' isochrone),',
             'if flag=2 there was a problem with spectroscopy,'
             'if flag=3 problem with photometry,'
             'if flag=4 problem with astrometry,'
             'if flag=5 problem with mass estimate.'
             ]
    return descr


def process_distances(input_data, output, id_col, name, which="Padova",
                      prior=True, npool=32,
                      with_parallax=True, with_mass=True):
    ''' Computes the distances to the stars in the table
        input_data. which is an array of strings which give the isochrone
        sets to use, if prior we use the Binney prior'''

    output_data = Table.from_pandas(input_data[[id_col]])
    output_data.meta['COMMENT'] = column_descr(name, id_col, which,
                                               prior, with_parallax, with_mass)

    output_data = new_distance_columns(output_data)
    output_data = new_add_columns(output_data)

    num_tasks = len(input_data)
    nindex = range(num_tasks)
    for Nindx in np.array_split(nindex, 100):
        if(npool > 1):
            p = Pool(npool)
            r = p.map(partial(process_single_distance, input_data, which,
                              prior, with_parallax, with_mass),
                      Nindx, chunksize=30)
            p.close()
            p.join()
        else:
            r = map(lambda i: process_single_distance(
                input_data, which, prior, with_parallax, with_mass, i), Nindx)
        for nindx, X in zip(Nindx, r):
            for K, n in enumerate(distance_output + additional_output):
                output_data[n][nindx] = X[K]
        output_data.write(output, path='data', format='hdf5',
                          compression=True,
                          serialize_meta=True,
                          overwrite=True)
    print 'All distances computed!'

    return output_data


def mfe(data):
    return np.nanmedian(data.e_fe_h)


def run_distance_pipeline(data, out, id_col, name, npool=-1):
    print mfe(data)
    isochrones = "Padova"
    init(types=isochrones,
         wemap=False, mean_feh_err=np.float64(mfe(data)))
    process_distances(data, out, id_col, name,
                      which=isochrones, prior=True,
                      npool=npool,
                      with_parallax=True,
                      with_mass=True)


if __name__ == '__main__':
    # run_distance_pipeline()
    init("Padova")
