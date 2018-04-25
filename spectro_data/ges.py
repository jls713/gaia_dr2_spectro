from utils import *
from spectro_data.utils import *
import sys
sys.path.append('../')
import cross_match


def format_columns(data):
    col_dict = {'japercormag4': 'Jv', 'hapercormag4': 'Hv',
                'kapercormag4': 'Kv',
                'japercormag4_err': 'eJv', 'hapercormag4_err': 'eHv',
                'kapercormag4_err': 'eKv',
                'RA': 'ra', 'DEC': 'dec',
                'TEFF': 'teff', 'E_TEFF': 'e_teff',
                'LOGG': 'logg', 'E_LOGG': 'e_logg',
                'FEH': 'fe_h', 'E_FEH': 'e_fe_h',
                'VRAD': 'hrv',
                'E_VRAD': 'e_hrv'}
    data = data.rename(index=str, columns=col_dict)
    fltr = (data['fe_h'] != data['fe_h'])
    data.loc[fltr, 'fe_h'] = data['MH'][fltr]
    data.loc[fltr, 'e_fe_h'] = data['E_MH'][fltr]
    return data


def load_data():
    ges = read_fits_to_pandas(
        '/data/jls/ges/GES_iDR4_WG15_Recommended_v2.fits')
    ges = ges[ges.RA == ges.RA].reset_index(drop=True)
    cm = crossmatch_2MASS(ges.RA.values, ges.DEC.values)
    for i in cm.keys():
        ges[i] = cm[i].values
    cv = crossmatch_vista(ges.RA.values, ges.DEC.values)
    for i in cv.keys():
        ges[i] = cv[i]

    ges = format_columns(ges)

    vista_fltr = (ges.Jv == ges.Jv)
    vista_h = (ges.Hv == ges.Hv)  # For VHS GPS

    ges['mag_use'] = [np.array(['Jv', 'Hv', 'Kv', 'G'])
                      for i in range(len(ges))]
    ges.loc[~vista_h, 'mag_use'] = \
        ges.loc[~vista_h].applymap(lambda x: np.array(['Jv', 'Kv', 'G']))
    ges.loc[~vista_fltr, 'mag_use'] = \
        ges.loc[~vista_fltr].applymap(lambda x: np.array(['J', 'H', 'K', 'G']))

    ges['rho_TZ'] = 0.
    ges['rho_gZ'] = 0.
    ges['rho_Tg'] = 0.

    ges['parallax'] = 0.
    ges['parallax_error'] = -1.
    ges['pmra'] = 0.
    ges['pmdec'] = 0.

    ges['mass'] = 0.
    ges['mass_error'] = -1.

    return ges


def load_and_match(output_file='/data/jls/GaiaDR2/spectro/GES_input.hdf5',
                   use_dr1=False):

    loaded, data = check_and_load(output_file, 'Gaia-ESO DR3')
    if loaded:
        return data

    data = load_data()
    data = cross_match.crossmatch_gaia_spectro(data, epoch=2000, dr1=use_dr1)

    write_input_file(data, output_file, 'Gaia-ESO DR3')

    return data
