###############################################################################
# APOGEE
##
# Most recent data: DR14
###############################################################################

# =============================================================================
# Load apogee-rave cross-match and find photometry (APASS, Gaia, WISE)
import numpy as np
import pandas as pd
from astropy.io import fits
# =============================================================================
from utils import *
sys.path.append('../')
import cross_match
sys.path.append('/data/jls/cyanide/comparisons/')
from neural_network import *
# =============================================================================

APOGEE_FOLDER = '/data/jls/apogee/apogee_data/'

DR12_file = 'allStar-v603-nodups.fits'
DR14_file = 'allStar-l31c.2.fits'


def load_data(calibrated=True, use_dr12=False):
    if use_dr12:
        fitsfile = APOGEE_FOLDER + DR12_file
    else:
        fitsfile = APOGEE_FOLDER + DR14_file

    apogeeF = fits.open(fitsfile)
    flds = ['APOGEE_ID',
            'J', 'H', 'K',
            'J_ERR', 'H_ERR', 'K_ERR',
            'RA', 'DEC',
            'ASPCAPFLAG', 'VHELIO_AVG', 'VSCATTER']
    apogeeR = apogeeF[1].data
    apogee = pd.DataFrame()
    for i in flds:
        apogee[i] = apogeeR[i].byteswap().newbyteorder()
    apogee['snr'] = apogeeR['SNR']
    apogee['AK'] = apogeeR['AK_TARG']
    apogee['lnAK'] = np.log(apogee['AK'])
    apogee.loc[apogee.AK <= 0., 'lnAK'] = np.log(
        np.min(apogee.AK[apogee.AK > 0.]))

    element_list = ['C', 'N', 'O', 'Mg', 'Al', 'Si', 'Ca', 'Fe', 'Ni', 'S']

    if not calibrated:
        param_fld = 'FPARAM'
        elem_fld = 'FELEM'
    else:
        param_fld = 'PARAM'
        elem_fld = 'X_H'
        if use_dr12:
            elem_fld = 'ELEM'

    params = ['TEFF', 'LOGG', 'M_H', 'ALPHA_M', 'C_M', 'N_M']
    params_indx = [0, 1, 3, 6, 4, 5]
    for i, f in zip(params_indx, params):
        apogee[f] = apogeeR[param_fld].T[i]
        apogee[f + '_ERR'] = np.sqrt(apogeeR[param_fld + '_COV'][:, i, i])
        apogee[f + '_FLAG'] = apogeeR['PARAMFLAG'].T[i]
    apogee['ASPCAPFLAG'] = apogeeR['ASPCAPFLAG']
    apogee['rho_TZ'] = apogeeR[param_fld + '_COV'][:, 0, 3] / \
        (apogee['M_H_ERR'] * apogee['TEFF_ERR'])
    apogee['rho_gZ'] = apogeeR[param_fld + '_COV'][:, 1, 3] / \
        (apogee['M_H_ERR'] * apogee['LOGG_ERR'])
    apogee['rho_Tg'] = apogeeR[param_fld + '_COV'][:, 0, 1] / \
        (apogee['TEFF_ERR'] * apogee['LOGG_ERR'])
    for i in ['rho_Tg', 'rho_gZ', 'rho_TZ']:
        fltr = (apogee[i] < -1.) | (apogee[i] > 1.)
        apogee.loc[fltr, i] = 0.

    if not use_dr12:
        aux_data = apogeeF[3].data
    else:
        aux_data = fits.open(fitsfile[:-12] + '.fits')[3].data
    elem_symbol = aux_data['ELEM_SYMBOL'][0]
    elem_index = [np.argwhere(elem_symbol == e)[0][0] for e in element_list]
    elem_symbol = elem_symbol[elem_index]
    # Is element stored relative to metallicity? -- note negation
    elem_met = 1 - aux_data['ELEMTOH'][0][elem_index]

    for i, e, em in zip(elem_index, element_list, elem_met):
        apogee[e + '_H'] = apogeeR[elem_fld].T[i] + \
            (not calibrated) * em * apogee['M_H']
        apogee[e + '_H_ERR'] = np.sqrt(apogeeR[elem_fld + '_ERR'][:, i]**2 +
                                       apogee['M_H_ERR'] * em *
                                       (not calibrated))
    elem_flg = (apogeeR['ELEMFLAG'][:, elem_index] > 0).T
    for e, ef in zip(elem_symbol, elem_flg):
        apogee[e + '_H_FLAG'] = ef

    apogee['C_N'] = apogee['C_H'] - apogee['N_H']
    apogee['C_N_ERR'] = np.sqrt(
        apogee['C_H_ERR']**2 + apogee['N_H_ERR']**2)

    if not use_dr12:
        # Remove duplicates
        apogee = apogee[apogeeR['EXTRATARG'] < 15].reset_index(drop=True)
        apogee = apogee.sort_values(
            'snr', ascending=False).reset_index(drop=True)
        apogee = apogee.drop_duplicates('APOGEE_ID')

    apogee['ID_2MASS'] = apogee['APOGEE_ID'].str.lstrip('2M')

    if use_dr12:
        for name in apogee.columns:
            if name not in flds and \
                    apogee[name].dtype.byteorder not in ('=', '|'):
                apogee[name] = apogee[name].values.byteswap().newbyteorder()

    apogee = apogee.replace(-9999., np.nan)

    apogee['mag_use'] = [np.array(['J', 'H', 'K', 'G'])
                            for i in range(len(apogee))]

    # Dwarfs -- assign broad logg = 4.5 \pm 1.5
    fltr = (apogee.TEFF == apogee.TEFF) & (
        apogee.M_H == apogee.M_H) & (apogee.LOGG != apogee.LOGG)
    apogee.loc[fltr, 'LOGG'] = 4.5
    apogee.loc[fltr, 'LOGG_ERR'] = 2.

    apogee['mass'] = 0.
    apogee['mass_error'] = -1.

    # Add masses
    with open('/data/jls/cyanide/comparisons/neural_network.pkl', 'r') as f:
        nn = pickle.load(f)
    flds = ['C_M', 'N_M', 'ALPHA_M', 'TEFF', 'LOGG', 'M_H']
    inputData = apogee[flds].values
    inputErrData = apogee[[fld + '_ERR' for fld in flds]].values
    results = nn.label(inputData, inputErrData, samples=50)
    fltr = True
    for f in flds:
        fltr &= (apogee[f] == apogee[f])
    apogee.loc[fltr, 'mass'], apogee.loc[fltr, 'mass_error'] = results[0][fltr], results[1][fltr]

    return apogee


def format_columns(data):
    col_dict = {'J_ERR': 'eJ', 'H_ERR': 'eH', 'K_ERR': 'eK',
                'M_H_ERR': 'e_fe_h', 'M_H': 'fe_h',
                'TEFF_ERR': 'e_teff', 'TEFF': 'teff',
                'LOGG_ERR': 'e_logg', 'LOGG': 'logg',
                'VHELIO_AVG': 'hrv',
                'VSCATTER': 'e_hrv',
                'RA': 'ra', 'DEC': 'dec'}
    data = data.rename(index=str, columns=col_dict)
    return data


def load_and_match(output_file='/data/jls/GaiaDR2/spectro/APOGEE_input.hdf5',
                   calibrated=True, use_dr12=False, use_dr1=False):

    loaded, data = check_and_load(output_file, 'APOGEE DR14')
    if loaded:
        return data

    data = load_data(calibrated=calibrated, use_dr12=use_dr12)
    data = format_columns(data)
    data = cross_match.crossmatch_gaia_spectro(data, dr1=use_dr1, epoch=2000.)

    write_input_file(data, output_file, 'APOGEE DR14')
    return data
