# compute extinction coefficients for passbands
try:
    from specutils.extinction import extinction_f99, ExtinctionF99, extinction_od94
except:
    pass
from astropy import units as u
from astropy.io import fits
import fitsio
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, splrep, splev
from scipy.integrate import quad
import glob
from schlafly2016 import *

##
# Schlafly et al. (2016) monotonic extinction curve
##


class schlafly_apogee(object):
    ''' Converts Green et al. (2017) extinction map values to A_lam but using
        Indebetouw (2005)calibration '''

    def __init__(self, RV=3.3):
        x = (RV - 3.3) / 9.1
        self.curve = extcurve(x)

    def __call__(self, lam):
        ''' angstrom '''
        return self.curve(lam) / 0.3118743701541266


class schlafly_apogee_green_calibrated(object):
    ''' Converts Green et al. (2017) extinction map values to A_lam '''

    def __init__(self, RV=3.3):
        x = (RV - 3.3) / 9.1
        self.curve = extcurve(x)

    def __call__(self, lam):
        ''' angstrom '''
        return (self.curve(lam) - 0.04121974) / 0.3118743701541266


def load_ck_model(teff=5000, logg=4.5, vega=False, schlafly=False):
    if vega:
        data = fits.open('/data/jls/stellar_models/ckm05_9500.fits')[1].data
        return data['WAVELENGTH'], data['g40']
    if schlafly:
        print('Not currently used but if used, double-check units of flux -- is it per lambda or per nu?')
        data = np.genfromtxt('/data/jls/stellar_models/munari_schlafly.dat')
        return data.T[0], data.T[1]
    data = fits.open('/data/jls/stellar_models/ckp00_' +
                     str(teff) + '.fits')[1].data
    s = '%0.1f' % logg
    return data['WAVELENGTH'], data['g%s%s' % (s[0], s[-1])]


filter_folder = ''  # '/data/jls/filters/'
filters = ['SDSSu', 'SDSSg', 'SDSSr', 'SDSSi', 'SDSSz',
           'LandoltB', 'LandoltV',
           'GaiaBP_dr2', 'GaiaG_dr2', 'GaiaRP_dr2', 'GaiaRVS',
           '2MASSJ', '2MASSH', '2MASSKs',
           'WISEW1', 'WISEW2',
           'PSg', 'PSr', 'PSi', 'PSz', 'PSy']
names = ['u', 'g', 'r', 'i', 'z',
         'B', 'V',
         'GBP', 'G', 'GRP', 'GRVS',
         'J', 'H', 'K',
         'W1', 'W2',
         'gP', 'rP', 'iP', 'zP', 'yP']
filter_dict = dict(zip(names, filters))


def load_filter(fltr):
    fltr_data = np.genfromtxt(filter_folder + filter_dict[fltr], skip_header=1)
    # Remove leading and trailing zeros
    fltr_data = fltr_data[fltr_data.T[1] > 0.]
    return fltr_data


def make_spectrum_figure():
    plt.figure(figsize=[17., 3.5])
    NN = 0
    ls = 'solid'

    # 1. Plot filters
    for fltr, n in zip(filters, names):
        fltr_data = load_filter(fltr)
        plt.plot(fltr_data.T[0] / 10.,  # nm
                 fltr_data.T[1] / np.max(fltr_data.T[1]), label=r'$%s$' % n,
                 ls=ls)
        NN += 1
        if(NN > 5):
            ls = 'dashed'
    plt.xlabel(r'$\lambda/\,\mathrm{nm}$')
    plt.ylabel('Response')
    plt.xlim(0., 5500.)
    plt.legend(loc='lower center', ncol=11, bbox_to_anchor=(0.5, 1.))

    # 2. Plot models
    ax2 = plt.gca().twinx()
    m = load_ck_model(teff=3500, logg=4.5)
    ax2.plot(m[0] / 10., np.log10(m[1] * 1e-2), color='k')
    m = load_ck_model(teff=6000, logg=4.5)
    ax2.plot(m[0] / 10., np.log10(m[1] * 1e-2), color='k')
    m = load_ck_model(teff=45000, logg=4.5)
    ax2.plot(m[0] / 10., np.log10(m[1] * 1e-2), color='k')
    m = load_ck_model(vega=True)
    ax2.plot(m[0] / 10., np.log10(m[1] * 1e-2), color='gray')
    ax2.set_xlim(0., 5500.)
    ax2.set_ylim(2., 10.)
    ax2.set_ylabel(r'Log[Flux/ $\mathrm{W\,m}^{-2}\mathrm{nm}^{-1}]$')
    plt.savefig('extinction_spectrum.pdf', bbox_inches='tight')


def lambda_eff(flt, m):
    f = load_filter(flt)
    f_interp = interp1d(f.T[0], f.T[1])
    m_interp = interp1d(m[0], m[1])
    upplim = np.max(f.T[0])
    if(upplim > 60000.):
        upplim = 59999.
    if(upplim > m[0][-1]):
        upplim = m[0][-1]
    return quad(lambda x: x**2 * m_interp(x) * f_interp(x),
                np.min(f.T[0]), upplim)[0]\
        / quad(lambda x: x * m_interp(x) * f_interp(x),
               np.min(f.T[0]), upplim)[0]

# ============================================================================
# Extinction coefficients Bayestar 17
# ============================================================================


def convolve_spectrum_band_schlafly(m, f, extinction_law,
                                    av=False, afactor=1.):

    m_interp = interp1d(m[0], m[1])
    f_interp = interp1d(f[0], f[1])

    def prod(x):
        aval = 1.
        if av:
            aval = np.power(10., -0.4 * extinction_law(x) * afactor)
        return f_interp(x) * m_interp(x) * aval * x
    upplim = np.max(f[0])
    if(upplim > 60000.):
        upplim = 59999.
    if(upplim > m[0][-1]):
        upplim = m[0][-1]
    return -2.5 * np.log10(quad(prod, np.min(f[0]), upplim)[0])

def compute_absolute_av_schlafly(m, RV=3.1, Afactor=1.):
    extinction_law = schlafly_apogee_green_calibrated(RV=RV)
    return extinction_law(lambda_eff('V',m))*Afactor

def compute_extinction_coeffs_schlafly(m, RV=3.3,
                                       filters_to_use=None,
                                       Afactor=0.01):
    extinction_law = schlafly_apogee_green_calibrated(RV=RV)
    results = {}
    use_filters = filter_dict
    if filters_to_use:
        use_filters = filters_to_use
    for f in use_filters:
        flt = load_filter(f)
        R = convolve_spectrum_band_schlafly(m, (flt.T[0], flt.T[1]),
                                            extinction_law, av=True,
                                            afactor=Afactor)
        R -= convolve_spectrum_band_schlafly(m, (flt.T[0], flt.T[1]),
                                             extinction_law, afactor=Afactor)
        results[f] = R / Afactor
    return results

# ============================================================================
# Extinction coefficients Bayestar 15
# ============================================================================


# It appears that the Fitzpatrick 1999 law is favoured
# (Schlafly 2010, 2011, Yuan 2013), although it appears
# to break down for UV (Yuan 2013).
# Schlafly et al. (2010) find $R_V=3.1\pm0.2$ and
# $N=0.78\pm0.06$.
# Here $N=A_{1\mu\mathrm{m,pred}}/(1.32 E(B-V)_\mathrm{SFD})$
# such that $E(B-V)=1.129 N E(B-V)_\mathrm{SFD} = 0.88062 E(B-V)_\mathrm{SFD}$
# i.e. Schlegel overpredicts by $\sim14\%$.


# AV = 1.  # Arbitrary choice -- needs to be small so in the linear regime
R_V, sig_R_V = 3.1, 0.2
N, sig_N = 0.78, 0.06
factor = R_V * 0.88062
try:
    EBV_Sch = extinction_od94(440.4 * u.nm, a_v=1., r_v=3.1) - \
    extinction_od94(542.8 * u.nm, a_v=1., r_v=3.1)
    Normratio = extinction_od94(1000. * u.nm, 1., r_v=3.1) / \
    extinction_f99(1000. * u.nm, 1., r_v=3.1) * N
    print('A_1micron/E(B-V)_SFD=', \
    extinction_od94(1000. * u.nm, a_v=1., r_v=3.1) / EBV_Sch)
except:
    pass

def convolve_spectrum_band(m, f, schlafly_factor, Extf99,
                           av=False, form='F99', afactor=1.):
    m_interp = interp1d(m[0], m[1])
    f_interp = interp1d(f[0], f[1])

    def prod(x):
        aval = 1.
        if av:
            aval = np.power(10., -0.4 * schlafly_factor * afactor *
                            Extf99((x / 10.) * u.nm)[()])
            # else:
            #     aval = np.power(
            #         10., -0.4 * extinction_od94((x / 10.) * u.nm,
            #                                     a_v=AV, r_v=3.1)[()])
        return f_interp(x) * m_interp(x) * aval * x
    upplim = np.max(f[0])
    if(upplim > 60000.):
        upplim = 59999.
    if(upplim > m[0][-1]):
        upplim = m[0][-1]
    return -2.5 * np.log10(quad(prod, np.min(f[0]), upplim)[0])

def compute_absolute_av(m, RV=3.1, Afactor=1.):
    Extf99 = ExtinctionF99(a_v=1., r_v=RV)
    schlafly_factor = extinction_od94(
        1000. * u.nm, a_v=1., r_v=3.1) / extinction_f99(1000. * u.nm,
                                                         a_v=1., r_v=RV) * N
    EBV = (extinction_od94(440.4 * u.nm, a_v=1., r_v=3.1) -
           extinction_od94(542.8 * u.nm, a_v=1., r_v=3.1))
    return schlafly_factor*Afactor*Extf99(lambda_eff('V',m)*0.1*u.nm)/EBV

def compute_extinction_coeffs(m, RV=3.1, filters_to_use=None, Afactor=0.01):
    ''' Afactor gives units of extinction curve to consider. Usually working
        in linear regime so we don't care about its value (provided it is
        small). To work out absolute values of monochrome extinction, we
        can use compute_absolute_av '''
    Extf99 = ExtinctionF99(a_v=1., r_v=RV)
    schlafly_factor = extinction_od94(
        1000. * u.nm, a_v=1., r_v=3.1) / extinction_f99(1000. * u.nm,
                                                         a_v=1., r_v=RV) * N
    EBV = (extinction_od94(440.4 * u.nm, a_v=1., r_v=3.1) -
           extinction_od94(542.8 * u.nm, a_v=1., r_v=3.1))
    schlafly_factor /= EBV

    results = {}
    use_filters = filter_dict
    if filters_to_use:
        use_filters = filters_to_use
    for f in use_filters:
        flt = load_filter(f)
        R = convolve_spectrum_band(m, (flt.T[0], flt.T[1]),
                                   schlafly_factor, Extf99, av=True,
                                   afactor=Afactor)\
            - convolve_spectrum_band(m, (flt.T[0], flt.T[1]),
                                     schlafly_factor, Extf99,
                                     afactor=Afactor)
        results[f] = R / Afactor # / EBV
    return results


# ============================================================================
# Extinction coefficients Bayestar 15
# ============================================================================


def write_extinction_table_latex(teff=4500, logg=4.5, RV=R_V,
                                 BSversion=2017):

    ext_file = 'extinction_coeff_%i_latex.dat' % BSversion
    if BSversion == 2017:
        extinction_coeff_fn = compute_extinction_coeffs_schlafly
    else:
        extinction_coeff_fn = compute_extinction_coeffs

    l_eff = {i: lambda_eff(i, load_ck_model(teff=teff, logg=logg))
             for i in filter_dict}
    results = extinction_coeff_fn(load_ck_model(teff=teff, logg=logg),
                                  RV=R_V)

    s = r'\begin{tabular}{lcccc}\\'
    s += r'Band&$\lambda_\mathrm{eff}/\,\AA$&$R(\lambda_\mathrm{eff})$&'
    s += r'$R(\lambda)$&$A_\lambda/A_V$&\\\hline '

    if BSversion == 2017:
        Data = {i: schlafly_apogee_green_calibrated(RV)(l_eff[i])
                for i in l_eff.keys()}
    else:
        Data = {i: extinction_f99(np.array(
            [l_eff[i] / 10.]) * u.nm, 1., r_v=RV)[0] * Normratio / EBV_Sch
            for i in l_eff.keys()}

    import operator
    for i in sorted(l_eff.items(), key=operator.itemgetter(1)):
        i = i[0]
        if i == 'G':
            pass
        s += r'%s&$%0.1f$&$%0.3f$&$%0.3f$&$%0.3f$\\' % (
            i, l_eff[i], Data[i], results[i], results[i]/results['V'])
    s += r'\hline\end{tabular}'
    with open(ext_file, 'w') as f:
        f.write(s)


def compute_extinction_table_lambdaeff(teff=4500, logg=4.5,
                                       BSversion=2017):

    ext_file = 'extinction_coeffs_%i_leff.dat' % BSversion
    l_eff = {i: lambda_eff(i, load_ck_model(teff=teff, logg=logg))
             for i in filter_dict}
    rvrange = np.linspace(2.1, 6.1, 21)
    if BSversion == 2017:
        Data = np.array([[schlafly_apogee_green_calibrated(R)(l_eff[i])
                          for R in rvrange]
                         for i in l_eff.keys()])
    else:
        Data = np.array([[extinction_f99(np.array(
            [l_eff[i] / 10.]) * u.nm, 1., r_v=R)[0] *
            Normratio / EBV_Sch for R in rvrange]
            for i in l_eff.keys()])
    with open(ext_file, 'w') as output:
        output.write(' '.join([str(r) for r in rvrange]) + '\n')
        [output.write(n + ' ' + ' '.join([str(dd) for dd in d]) + '\n')
         for n, d in zip(l_eff.keys(), Data)]


def compute_extinction_table(teff=4500, logg=4.5,
                             BSversion=2017):
    ext_file = 'extinction_coeffs_%i.dat' % BSversion

    if BSversion == 2017:
        extinction_coeff_fn = compute_extinction_coeffs_schlafly
    else:
        extinction_coeff_fn = compute_extinction_coeffs

    rvrange = np.linspace(2.1, 6.1, 21)
    Data = [extinction_coeff_fn(load_ck_model(
        teff=teff, logg=logg), RV=rv) for rv in rvrange]
    pData = np.array([[Data[rv][m] for rv in range(len(rvrange))]
                      for m in filter_dict.keys()])
    with open(ext_file, 'w') as output:
        output.write(' '.join([str(r) for r in rvrange]) + '\n')
        [output.write(n + ' ' + ' '.join([str(dd) for dd in d]) + '\n')
         for n, d in zip(filter_dict.keys(), pData)]


def compute_Gaia_extinction_grid(BSversion=2017):
    ext_file = 'extinction_coeffs_G_%s.dat' % BSversion

    if BSversion == 2017:
        extinction_coeff_fn = compute_extinction_coeffs_schlafly
    else:
        extinction_coeff_fn = compute_extinction_coeffs

    files = glob.glob('/data/jls/stellar_models/ckp00*')
    rvrange = np.linspace(2.1, 6.1, 21)
    ResultsG = np.zeros((len(files), len(rvrange)))
    TEff = np.zeros(len(files))
    for n, f in enumerate(files):
        teff = f.split('_')[-1][:-5]
        TEff[n] = teff
        for nj, j in enumerate(rvrange):
            ResultsG[n][nj] = extinction_coeff_fn(load_ck_model(
                teff=teff, logg=4.5), RV=j, filters_to_use=['G'])['G']
        print(teff)
    with open(ext_file, 'w') as output:
        output.write(' '.join([str(np.log10(t))
                               for t in TEff[np.argsort(TEff)][:-1]]) + '\n')
        [output.write(' '.join([str(dd) for dd in d]) + '\n')
         for d in ResultsG[np.argsort(TEff)][:-1]]


def load_extinction_data(version=2017):
    redd_maps = {}
    with open('extinction_coeffs_%i.dat' % version) as f:
        R_V_grid = np.fromstring(f.readline(), dtype=np.float64, sep=' ')
        for l in f.readlines():
            spl = l.split(" ", 1)
            redd_maps[spl[0]] = np.fromstring(
                spl[1], dtype=np.float64, sep=' ')

    interp_maps = {n: interp1d(R_V_grid, redd_maps[n]) for n in redd_maps}

    R_G = np.zeros((75, 21))
    with open('extinction_coeffs_G_%i.dat' % version) as f:
        logTeffgrid = np.fromstring(f.readline(), dtype=np.float64, sep=' ')
        for n, l in enumerate(f.readlines()):
            R_G[n] = np.fromstring(l, dtype=np.float64, sep=' ')

    from scipy.interpolate import interp2d
    interp_G_maps = interp2d(R_V_grid, logTeffgrid, R_G)
    return R_V_grid, logTeffgrid, R_G, interp_maps, interp_G_maps

# ============================================================================
# G, GRP and GBP from models
# ============================================================================


def compute_gaiamags(model, mode='dr2'):
    vega = load_ck_model(vega=True)
    fls = {'dr2': 'GaiaDR2_ZeroPoints.dat',
           'revised': 'GaiaDR2_RevisedZeroPoints.dat',
           'premission': 'na.dat'}
    zp = np.genfromtxt(fls[mode])
    zp = {'G': zp[0][0], 'GBP': zp[0][2], 'GRP': zp[0][4]}
    mag = {}
    for fltr in ['G', 'GBP', 'GRP']:
        flt = load_filter(fltr)
        r = convolve_spectrum_band_schlafly(model, (flt.T[0], flt.T[1]),
                                            None, av=False, afactor=0.)
        r -= convolve_spectrum_band_schlafly(vega, (flt.T[0], flt.T[1]),
                                             None, av=False, afactor=0.)
        r += zp[fltr]
        mag[fltr] = r
    return mag


def gaiamag_with_effT(mode='dr2'):
    files = glob.glob('/data/jls/stellar_models/ckp00*')
    teffrange = np.zeros(len(files), dtype=int)
    G = np.zeros(len(teffrange))
    GBP = np.zeros_like(G)
    GRP = np.zeros_like(G)
    for i, f in enumerate(files):
        teffrange[i] = f.split('_')[-1][:-5]
        model = load_ck_model(teff=teffrange[i])
        mag = compute_gaiamags(model)
        G[i], GBP[i], GBP[i] = mag['G'], mag['GBP'], mag['GRP']
    args = np.argsort(teffrange)
    teffrange = teffrange[args]
    G = G[args]
    GBP = GBP[args]
    GRP = GRP[args]
    return teffrange, G, GBP, GRP

# ============================================================================
# Non-linearity
# ============================================================================


def find_aiav_gradient(model, fltr, BSversion=2017, RV=3.1):
    if BSversion == 2017:
        extinction_coeff_fn = compute_extinction_coeffs_schlafly
        compute_absolute = compute_absolute_av_schlafly
    else:
        extinction_coeff_fn = compute_extinction_coeffs
        compute_absolute = compute_absolute_av

    Afactor = np.logspace(-2., 0., 4)
    ff = [extinction_coeff_fn(model, filters_to_use=[fltr],
                              Afactor=af, RV=RV)[fltr] for af in Afactor]
    AA = np.array([extinction_coeff_fn(model, filters_to_use=['V'],
                              Afactor=af, RV=RV)['V'] for af in Afactor])*Afactor
    # AA = [compute_absolute(model, RV=RV, Afactor=af) for af in Afactor]
    m, b = np.polyfit(AA, ff / ff[0], 1)

    return m
# ============================================================================
# Polynomial fits
# ============================================================================


def poly2d(coeff, x, y):
    return coeff[0] + coeff[1] * x + coeff[2] * y + coeff[3] * x * x +\
        coeff[4] * x * x * y + coeff[5] * y * y + coeff[6] * x * y * y +\
        coeff[7] * x * y +\
        coeff[8] * x * x * x + coeff[9] * y * y * y


def polynomial_fit_G_extinction(fltr_set='2MASS', version=2017, RV=3.1,
                                block=False):

    rvrange, logTeffgrid, R_G, interp_maps, interp_G_maps \
        = load_extinction_data(version)

    # Using contourf to provide my colorbar info, then clearing the figure
    Z = [[0, 0], [0, 0]]
    CS3 = plt.contourf(Z,
                       np.append(
                           rvrange - 0.5 * (rvrange[1] - rvrange[0]),
                           rvrange[-1] + 0.5 * (rvrange[1] - rvrange[0])),
                       cmap=plt.cm.plasma)
    plt.clf()
    plt.figure(figsize=[3.32,2.5])

    iso_fl = {'2MASS': '2mass_wise_solar.dat',
              'SDSS': 'sdss_solar.dat',
              'Landolt': 'ubvrijhk_solar.dat',
              'Pan-STARRS': 'panstarrs_solar.dat',
              'GBP_GRP': 'gaia_solar.dat',
              'G_GRP': 'gaia_solar.dat',
              'GBP_G': 'gaia_solar.dat'}
    iso_in = {'2MASS': [8, 10], 'SDSS': [9, 10], 'Landolt': [9, 10],
              'Pan-STARRS': [8, 9],
              'GBP_GRP': [9, 10],
              'G_GRP': [8, 10],
              'GBP_G': [9, 8]}
    iso_lbl = {'2MASS': r'$(J-K_s)$', 'SDSS': r'$(g-r)$',
               'Landolt': r'$(B-V)$',
               'Pan-STARRS': r'$(g_P-r_P)$',
               'GBP_GRP': r'$(G_{BP}-G_{RP})$',
               'G_GRP': r'$(G-G_{RP})$',
               'GBP_G': r'$(G_{BP}-G)$',
               'Teff': r'$\log_{10}(T_\mathrm{eff})-4$'}

    if fltr_set in iso_fl.keys():
        tmass_next = np.genfromtxt(
            '/data/jls/isochrones/ob_stars/isochrones/' + iso_fl[fltr_set])
        x = tmass_next.T[iso_in[fltr_set][0]] - \
            tmass_next.T[iso_in[fltr_set][1]]
        y = tmass_next.T[5]

        fltr = True#(tmass_next.T[6]<4.5)
        print('WITH LOGG CUTS')
        x, y = x[fltr], y[fltr]
    elif fltr_set is 'Teff':
        x = logTeffgrid - 4.
        y = logTeffgrid - 4.
    else:
        teff, g, gbp, grp = gaiamag_with_effT()
        y = np.log10(teff)
        fltr_set_spl = fltr_set.split('_')
        ddict = {'G': g, 'GBP': gbp, 'GRP': grp}
        x = ddict[fltr_set_spl[0]] - ddict[fltr_set_spl[1]]

    if fltr_set is 'Teff':
        poly_f = lambda x: x - 4.
    else:
        z = np.polyfit(y[(y < 5.) & (y > 3.47)], x[(y < 5.) & (y > 3.47)], 8)
        poly_f = np.poly1d(z)
        x = x[(y < 5.) & (y > 3.47)]

    # print logTeffgrid, np.power(10.,logTeffgrid)

    plt.xlim(np.min(x), np.max(x))
    plt.ylim(1.5, 3.5)
    for rv in rvrange:
        plt.plot(poly_f(logTeffgrid), interp_G_maps(
            rv, logTeffgrid),
            color=plt.cm.plasma((rv - np.min(rvrange)) /
                                (np.max(rvrange) - np.min(rvrange))))
    plt.colorbar(CS3, label=r'$R_V$')
    plt.xlabel(iso_lbl[fltr_set])
    plt.ylabel(r'$R(G)$')
    plt.axhline(2.977, color='k', ls='dashed')
    plt.axhline(2.508, color='k', ls='dashed')
    plt.annotate(r'$R(r)$', xy=(1., 2.6), fontsize=14)
    plt.annotate(r'$R(V)$', xy=(1., 3.1), fontsize=14)

    X, Y = np.meshgrid(poly_f(logTeffgrid), rvrange - RV, indexing='ij')
    Z = R_G

    X = X.flatten()
    Y = Y.flatten()
    A = np.array([X * 0 + 1, X, Y, X**2, X**2 * Y, Y **
                  2, X * Y**2, X * Y, X**3, Y**3]).T
    B = Z.flatten()

    coeff, r, rank, s = np.linalg.lstsq(A, B)
    z = np.array([poly2d(coeff, 0., yy) for xx, yy in zip(X, Y)])
    # plt.plot(Y,z,'.',ms=8,color='k',label='Polynomial fit, $R_V=3.1$')

    poly1d = np.polyfit(poly_f(logTeffgrid),
                        interp_G_maps(RV, logTeffgrid), 3).flatten()
    plt.plot(poly_f(logTeffgrid)[::-1],
             np.poly1d(poly1d)(poly_f(logTeffgrid)[::-1]), '.-', ms=4,
             color='k', label='Polynomial fit, $R_V=%0.1f$' % RV)
    plt.legend(loc='lower center', bbox_to_anchor=(0.5, 0.9))
    print(fltr_set, ','.join('%0.4f' % f for f in coeff))
    print(fltr_set, ','.join('%0.4f' % f for f in poly1d[::-1]))

    with open(fltr_set + '_' + str(version) + '.poly', 'w') as fll:
        fll.write(fltr_set + '\n')
        fll.write(','.join('%0.4f' % f for f in coeff) + '\n')
        fll.write(','.join('%0.4f' % f for f in poly1d[::-1]) + '\n')
    plt.show(block=block)
    #plt.savefig('extinction_Gaia_JK.pdf', bbox_inches='tight')


if __name__ == '__main__':

    # Bayestar 2017
    RV = 3.3
    #write_extinction_table_latex(RV=RV)
    #compute_extinction_table_lambdaeff()
    #compute_extinction_table()
    #exit()
    #compute_Gaia_extinction_grid()
    #polynomial_fit_G_extinction(fltr_set='Teff', RV=RV)
    #polynomial_fit_G_extinction(fltr_set='2MASS', RV=RV, )
    #polynomial_fit_G_extinction(fltr_set='SDSS', RV=RV)
    #polynomial_fit_G_extinction(fltr_set='Landolt', RV=RV)
    #polynomial_fit_G_extinction(fltr_set='Pan-STARRS', RV=RV)
    #polynomial_fit_G_extinction(fltr_set='G_GRP', RV=RV)
    #polynomial_fit_G_extinction(fltr_set='GBP_G', RV=RV)
    #polynomial_fit_G_extinction(fltr_set='GBP_GRP', RV=RV)
    #model = load_ck_model(teff=4500)
    #t4500 = [find_aiav_gradient(model, 'G', BSversion=2017, RV=RV),
    #         find_aiav_gradient(model, 'GBP', BSversion=2017, RV=RV),
    #         find_aiav_gradient(model, 'GRP', BSversion=2017, RV=RV)]
    #model = load_ck_model(teff=10000)
    #t10000 = [find_aiav_gradient(model, 'G', BSversion=2017, RV=RV),
    #          find_aiav_gradient(model, 'GBP', BSversion=2017, RV=RV),
    #          find_aiav_gradient(model, 'GRP', BSversion=2017, RV=RV)]
    #np.savetxt('BS2017_nonlinearity.dat', np.vstack((t4500, t10000)))

    # # Bayestar 2015
    RV = 3.1
    write_extinction_table_latex(BSversion=2015, RV=RV)
    compute_extinction_table_lambdaeff(BSversion=2015)
    compute_extinction_table(BSversion=2015)
    # compute_Gaia_extinction_grid(BSversion=2015)
    # polynomial_fit_G_extinction(fltr_set='Teff', RV=RV, version=2015)
    # polynomial_fit_G_extinction(fltr_set='2MASS', RV=RV, version=2015)
    # polynomial_fit_G_extinction(fltr_set='SDSS', RV=RV, version=2015)
    # polynomial_fit_G_extinction(fltr_set='Landolt', RV=RV, version=2015)
    # polynomial_fit_G_extinction(fltr_set='Pan-STARRS', RV=RV, version=2015)
    # polynomial_fit_G_extinction(fltr_set='G_GRP', RV=RV, version=2015)
    # polynomial_fit_G_extinction(fltr_set='GBP_G', RV=RV, version=2015)
    # polynomial_fit_G_extinction(fltr_set='GBP_GRP', RV=RV, version=2015)
    # model = load_ck_model(teff=4500)
    # t4500 = [find_aiav_gradient(model, 'G', BSversion=2015),
    #          find_aiav_gradient(model, 'GBP', BSversion=2015),
    #          find_aiav_gradient(model, 'GRP', BSversion=2015)]
    # model = load_ck_model(teff=10000)
    # t10000 = [find_aiav_gradient(model, 'G', BSversion=2015),
    #           find_aiav_gradient(model, 'GBP', BSversion=2015),
    #           find_aiav_gradient(model, 'GRP', BSversion=2015)]
    # np.savetxt('BS2015_nonlinearity.dat', np.vstack((t4500, t10000)))
