import numpy as np
import isodist_js as isodist

which_iso = "BaSTI"

# Initialize a grid of BaSTI VISTA band isochrones and load extinction map


def init():
    isodist.init_isochrone(which_iso, 1, 0.4)
    isodist.load_emap()
    isodist.load_prior("2018")


if __name__ == '__main__':

    init()

    # Fake data
    mag_str = np.array(["Hv", "Jv"])
    mags = np.array([16.6463, 16.8812])
    Z, TEFF, LOGG = -0.35, 3.7751, 4.55148
    ERRZ, ERRTEFF, ERRLOGG = 0.02, 0.005, 0.06
    l, b = 0.4, 1.  # in radians
    errmags = np.array([0.022, 0.022])

    # Get data in correct format
    data = np.array([np.float64(Z), np.float64(TEFF),
                     np.float64(LOGG), np.float64(l), np.float64(b)])
    data_errs = np.array(
        [np.float64(ERRZ), np.float64(ERRTEFF), np.float64(ERRLOGG)])

    prior = True  # Use Binney spatial prior
    print isodist.prob_distance(mags, data, errmags, data_errs,
                                prior, which_iso, mag_str, True, 0.,
                                False, 0., -1.)
