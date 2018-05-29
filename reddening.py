import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
import astropy.units as units
from dustmaps.bayestar import BayestarQuery
from dustmaps.marshall import MarshallQuery
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import quad
from scipy.interpolate import RegularGridInterpolator
import mwdust
import healpy
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))+'/'
# with open('setup.json') as data_file:
#     data = json.load(data_file)
#     isochrone_folder=data['isochrone_folder']

# All E(B-V) in Bayestar 2017 units

AVconst_BS2015 = 2.742
AVconst_BS2018 = 2.9492252258418326

class ReddeningMaps(object):
    def __init__(self, with_extinction_maps=False):

        self.redd_maps = {}
        with open(dir_path+'extinction/extinction_coeffs_2017.dat') as f:
            self.R_V_grid = np.fromstring(
                f.readline(), dtype=np.float64, sep=' ')
            for l in f.readlines():
                spl = l.split(" ", 1)
                self.redd_maps[spl[0]] = np.fromstring(
                    spl[1], dtype=np.float64, sep=' ')

        self.interp_maps = {n: interp1d(
            self.R_V_grid, self.redd_maps[n]) for n in self.redd_maps}

        self.R_G = np.zeros((75, len(self.R_V_grid)))
        with open(dir_path+'extinction/extinction_coeffs_G_2017.dat') as f:
            self.logTeffgrid = np.fromstring(
                f.readline(), dtype=np.float64, sep=' ')
            for n, l in enumerate(f.readlines()):
                self.R_G[n] = np.fromstring(l, dtype=np.float64, sep=' ')

        self.interp_G_maps = interp2d(
            self.R_V_grid, self.logTeffgrid, self.R_G)

        if(with_extinction_maps):
            self.sfd = SFDQuery()
            self.bayestar = BayestarQuery(max_samples=10)
            self.marshall = MarshallQuery()
            self.drimmel = mwdust.Drimmel03(filter='Landolt V')

    def get_log_ebv_along_los(self, ldeg, bdeg, RV=3.3):

        median_log_ebv, std_log_ebv, dist_ebv = \
            self.get_bayestar_log_ebv_along_los(ldeg, bdeg)
        # Invalid value encountered here means it is outside the footprint
        in_bayestar_footprint = (median_log_ebv[0] == median_log_ebv[0])
        if not in_bayestar_footprint:
            # Use Marshall map (2006) -- returns A_Ks
            distance_grid = np.logspace(-1.2, 1.8, 31) * units.kpc
            median_log_ebv, std_log_ebv, dist_ebv = \
                self.get_marshall_log_ebv_along_los(ldeg, bdeg, distance_grid,
                                                    RV=RV)
            in_marshall_footprint = (median_log_ebv[0] == median_log_ebv[0])
            if not in_marshall_footprint:
                median_log_ebv, std_log_ebv, dist_ebv = \
                    self.get_drimmel_log_ebv_along_los(ldeg, bdeg,
                                                       distance_grid,
                                                       RV=RV)
        return median_log_ebv.astype(np.float32), \
            std_log_ebv.astype(np.float32), \
            dist_ebv.astype(np.float32)

    def get_log_av_along_los(self, l, b, RV=3.3):
        median_log_ebv, std_log_ebv, dist_ebv = \
            self.get_log_ebv_along_los(l, b)
        AV = self.get_reddening('V', RV=RV)
        median_log_ebv += np.log(AV)
        return median_log_ebv, std_log_ebv, dist_ebv

    def write_extinction_map(self, nside=1024, RV=3.3):
        def pix2lb(pix, nest=False):
            """Convert pixel ID to Galactic l and b"""
            __, _ = healpy.pix2ang(nside, pix, nest=nest)
            l, b = np.rad2deg(_), np.rad2deg(np.pi / 2 - __)
            return l, b
        pix_ids = np.arange(healpy.nside2npix(nside))
        lb = pix2lb(pix_ids)
        distance_grid = np.logspace(-1.2, 1.8, 31)
        coords = [SkyCoord(lb[0] * units.deg, lb[1] * units.deg,
                           distance=d * units.kpc,
                           frame='galactic') for d in distance_grid]
        bayestar = np.array([np.log(self.bayestar(c, mode='median'))
                             for c in coords]).T
        marshall = np.array([np.log(self.marshall(c)) for c in coords]).T
        marshall -= np.log(self.get_reddening('K', RV=RV))
        drimmel = np.log(self.drimmel.evaluate_vector(lb[0], lb[1],
                                                      distance_grid))
        drimmel -= np.log(self.get_reddening('V', RV=RV))
        bayestar[bayestar != bayestar] = marshall[bayestar != bayestar]
        bayestar[bayestar != bayestar] = drimmel[bayestar != bayestar]
        # return bayestar
        with open("config.json") as config:
            config = json.load(config)
        healpy.fitsfunc.write_map(
            config["dir"]["extinction"] + "combined_nside%i.fits" % nside,
            bayestar.T, nest=False, overwrite=True)

    def get_bayestar_ebv(self, l, b, s):
        coords = SkyCoord(l * units.deg, b * units.deg,
                          distance=s * units.kpc, frame='galactic')
        ebv_samples = self.bayestar(coords, mode='samples')
        with np.errstate(invalid='ignore'):
            return np.median(ebv_samples, axis=-1), \
                np.std(ebv_samples, axis=-1)

    def get_bayestar_ebv_along_los(self, l, b):
        coords = SkyCoord(l * units.deg, b * units.deg, frame='galactic')
        ebv_samples = self.bayestar(coords, mode='samples')
        with np.errstate(invalid='ignore'):
            return np.median(ebv_samples, axis=0), \
                np.std(ebv_samples, axis=0), self.bayestar.distances

    def get_bayestar_log_ebv_along_los(self, l, b):
        coords = SkyCoord(l * units.deg, b * units.deg, frame='galactic')
        ebv_samples = self.bayestar(coords, mode='samples')
        with np.errstate(invalid='ignore'):
            return np.median(np.log(ebv_samples), axis=0), \
                np.std(np.log(ebv_samples), axis=0), self.bayestar.distances

    def get_marshall_log_ebv_along_los(self, l, b, s, RV=3.3):
        coords = SkyCoord(l * units.deg, b * units.deg,
                          distance=s, frame='galactic')
        ak = self.marshall(coords, return_sigma=True)
        ebv = ak[0] / self.get_reddening('K', RV=RV)  # EBV
        return np.log(ebv), ak[1] / ak[0], s

    def get_drimmel_log_ebv_along_los(self, l, b, s, RV=3.3,
                                      assumed_sigma=0.3):
        '''
            Currently doesnt work for vector l,b  but in mwdust.Drimmel03 there
            is available code
        '''
        av = self.drimmel(l, b, s.value)
        ebv = av / self.get_reddening('V', RV=RV)  # EBV
        return np.log(ebv), assumed_sigma * np.ones_like(ebv), s

    def get_ebv(self, l, b):
        ''' In Bayestar 2017 units '''
        ebv = self.sfd(SkyCoord(l, b, unit='deg', frame='galactic'))
        return ebv * AVconst_BS2015 / AVconst_BS2018

    def plot_schlegel(self, ax, llims=[0., 360.], blims=[-90., 90.]):
        ''' Plot the Schlegel E(B-V) map '''
        ll = np.linspace(llims[0], llims[1], 300)
        b = np.linspace(blims[0], blims[1], 300)
        ll, b = np.meshgrid(ll, b)
        coords = SkyCoord(ll, b, unit='deg', frame='galactic')
        AV = self.sfd(coords) * AVconst_BS2015 / AVconst_BS2018
        ax.imshow(np.log(AV),
                  origin='lower',
                  interpolation='nearest',
                  extent=[llims[0], llims[1], blims[0], blims[1]],
                  cmap=plt.cm.binary)
        ax.set_xlim(0., 360.)
        ax.set_ylim(-90., 90.)

    def get_reddening(self, c, RV=3.3, logTeff=4.):
        if c == 'G':
            return self.interp_G_maps(RV, logTeff)
        else:
            return self.interp_maps[c](RV)

    def get_reddening_color(self, c, RV=3.3, logTeff=4.):
        return self.get_reddening(c[0], RV, logTeff) - \
            self.get_reddening(c[1], RV, logTeff)

    def reddening_vector(self, color1, color2, RV=3.3, logTeff=4.):
        return [self.get_reddening_color(color1, RV, logTeff),
                self.get_reddening_color(color2, RV, logTeff)]

    def color_along_reddening_line(self, color, color1, color2,
                                   offy=0., offx=0., RV=3.3, logTeff=4.):
        rr = self.reddening_vector(color1, color2, RV=RV, logTeff=logTeff)
        return (color - offx) * rr[1] / rr[0] + offy

    def plot_reddening_vector(self, color1, color2, label="",
                              offy=0., offx=0., RV=3.3, logTeff=4.):
        xx = np.linspace(0., 7.)
        rr = self.reddening_vector(color1, color2, RV=RV, logTeff=logTeff)
        return plt.plot(xx * rr[0] + offx, xx * rr[1] + offy,
                        color='r', label=label, lw=3)

    def ak_to_ebv(self, AK, RV=3.3):
        return AK / self.get_reddening('K', RV)


def integrate_toy_extinction_map(l, b, s, ebv):
    def density(s):
        R0 = 8.3
        hR = 4.2
        hz = 0.088
        gammafl = 0.0054
        Rw = 8.4
        Rfl = 1.12 * R0
        gammaw = 0.18
        R = np.sqrt((R0 - s * np.cos(l) * np.cos(b)) **
                    2 + (s * np.sin(l) * np.cos(b))**2)
        sinp = s * np.sin(l) * np.cos(b) / R
        z = s * np.sin(b)
        kfl = 1. + gammafl * min(Rfl, R - Rfl)
        zw = gammaw * min(Rw, R - Rw) * sinp
        return np.exp((R0 - R) / hR - np.fabs(z - zw) / kfl / hz)
    epsabs = 1e-3
    epsrel = 1e-3
    return np.array([ebv * quad(density, 0, ss,
                                epsabs=epsabs, epsrel=epsrel)[0] for ss in s])\
        / quad(density, 0, np.inf, epsabs=epsabs, epsrel=epsrel)[0]


class combo_extinction_map(ReddeningMaps):
    def __init__(self):
        ReddeningMaps.__init__(self, with_extinction_maps=True)

    def __call__(self, l, b, s):
        median_ebv, std_ebv = self.get_bayestar_ebv(l, b, s)
        not_in_bayestar_footprint = (median_ebv != median_ebv)
        median_ebv[not_in_bayestar_footprint] = map(
            lambda L: integrate_toy_extinction_map(
                np.deg2rad(L[0]),
                np.deg2rad(L[1]),
                s, self.get_ebv(L[0], L[1])),
            zip(l[not_in_bayestar_footprint],
                b[not_in_bayestar_footprint]))
        std_ebv[not_in_bayestar_footprint] = \
            median_ebv[not_in_bayestar_footprint] / 2.
        return median_ebv, std_ebv


class gridded_combo_extinction_map(ReddeningMaps):
    def __init__(self, NL=60, NB=30):
        ReddeningMaps.__init__(self, with_extinction_maps=True)
        sgrid = np.log(
            self.get_bayestar_log_ebv_along_los(np.array([0.]),
                                                np.array([0.]))[-1].value)

        def ext(l, b, s):
            median_ebv, std_ebv = self.get_bayestar_ebv(l, b, s)
            not_in_bayestar_footprint = (median_ebv != median_ebv)
            median_ebv[not_in_bayestar_footprint] = map(
                lambda L: integrate_toy_extinction_map(
                    np.deg2rad(L[0]),
                    np.deg2rad(L[1]),
                    s, self.get_ebv(L[0], L[1])),
                zip(l[not_in_bayestar_footprint],
                    b[not_in_bayestar_footprint]))
            std_ebv[not_in_bayestar_footprint] = \
                median_ebv[not_in_bayestar_footprint] / 2.
            return median_ebv, std_ebv
        lgrid = np.linspace(0., 360., NL)
        bgrid = np.linspace(-90., 90., NB)
        L, B, S = np.meshgrid(lgrid, bgrid, sgrid, indexing='ij')

        ff = np.array([map(lambda l:ext(np.array([l[0]]),
                                        np.array([l[1]]),
                                        np.exp(np.array([l[2]]))),
                           zip(L.flatten(), B.flatten(), S.flatten()))])
        self.interp_med = RegularGridInterpolator((lgrid, bgrid, sgrid),
                                                  np.reshape(
                                                      ff[0, :, 0, 0],
                                                      np.shape(L)),
                                                  fill_value=None)
        self.interp_std = RegularGridInterpolator((lgrid, bgrid, sgrid),
                                                  np.reshape(
                                                      ff[0, :, 1, 0],
                                                      np.shape(L)),
                                                  fill_value=None)

    def __call__(self, l, b, s):
        lbs = np.vstack((l, b, np.log(s))).T
        return self.interp_med(lbs), self.interp_std(lbs)


if __name__ == '__main__':
    # print integrate_toy_extinction_map(0.3, 0.3, np.array([1., 2.]), 1.)
    red = ReddeningMaps(with_extinction_maps=True)
    red.write_extinction_map()
    # print red.get_bayestar_ebv_along_los(0., 0.)
