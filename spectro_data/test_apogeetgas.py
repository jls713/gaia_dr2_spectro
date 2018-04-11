###############################################################################
# APOGEE-TGAS test set from Payel
###############################################################################

# =============================================================================
# Load apogee-tgas test set
import numpy as np
import pandas as pd
# =============================================================================
from utils import *
import sys
sys.path.append('/data/jls/cyanide/comparisons/')
sys.path.append('../')
import cross_match
from neural_network import *
# =============================================================================

APOGEE_FOLDER = '/data/jls/GaiaDR2/spectro/'

data_file = 'APOGEE_TGAS_DR14_supp_keplercannon_masses_ages.csv'


def variance_weighted_alpha(labelled_set):
    cols = ['[Mg/M]', '[Si/M]', '[O/M]', '[Ca/M]']
    weights = np.array([1. / np.median(labelled_set['err_' + s][
                                       (labelled_set['err_' + s] > 0.)]**2)
                        for i, s in enumerate(cols)])
    err = np.sqrt(1. / np.sum(weights))
    weights /= np.sum(weights)
    return (labelled_set['[Mg/M]'] * weights[0] +
            labelled_set['[Si/M]'] * weights[1] +
            labelled_set['[O/M]'] * weights[2] +
            labelled_set['[Ca/M]'] * weights[3]), err


def load_data():

    apogee = pd.read_csv(APOGEE_FOLDER + data_file)
    apogee = apogee.replace(-9999., np.nan)

    apogee['mag_use'] = [np.array(['J', 'H', 'K'])
                         for i in range(len(apogee))]

    apogee['mass'] = apogee['kepmass']
    apogee['mass_error'] = .5 * (apogee['kepmass68U'] - apogee['kepmass68L'])

    # No correlations
    apogee['rho_TZ'] = 0.
    apogee['rho_Tg'] = 0.
    apogee['rho_gZ'] = 0.

    # Add masses
    with open('/data/jls/cyanide/comparisons/neural_network.pkl', 'r') as f:
        nn = pickle.load(f)
    apogee['ALPHA_M'], apogee['err_ALPHA_M'] = variance_weighted_alpha(apogee)
    inputData = apogee[['[C/M]', '[N/M]', 'ALPHA_M',
                        'teff', 'logg', 'm_h']].values
    inputErrData = apogee[['err_[C/M]', 'err_[N/M]', 'err_ALPHA_M',
                           'teff_err', 'logg_err', 'm_h_err']].values
    apogee['mass'], apogee['mass_error'] = nn.label(inputData, inputErrData)

    return apogee


def format_columns(data):
    col_dict = {'j': 'J', 'h': 'H', 'k': 'K',
                'ej': 'eJ', 'eh': 'eH', 'ek': 'eK',
                'm_h_err': 'e_fe_h', 'm_h': 'fe_h',
                'teff_err': 'e_teff', 'teff': 'teff',
                'logg_err': 'e_logg', 'logg': 'logg',
                'vlos': 'hrv'}
    data = data.rename(index=str, columns=col_dict)
    return data


def load_and_match(output_file='/data/jls/GaiaDR2/spectro/APOGEE_TGAS_test.hdf5'):

    loaded, data = check_and_load(output_file, 'APOGEE TGAS test set')
    if loaded:
        return data

    data = load_data()
    data = format_columns(data)
    data = cross_match.crossmatch_gaia_spectro(data, dr1=True, epoch=2000.)

    write_input_file(data, output_file, 'APOGEE TGAS test set')
    return data
