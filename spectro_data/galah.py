from spectro_data.utils import *
sys.path.append('../')
import cross_match


def load_data_dr1():
    data = pd.read_csv('/data/jls/galah/galah.dat', sep=r'\s+',
                       skipinitialspace=True,
                       names=['GALAH_ID', 'TYCHO2ID', '2MASSID',
                              'ra', 'dec', 'teff', 'logg', 'fe_h',
                              'alpha_fe', 'hrv', '(m-M)v', 'E(B-V)'])
    data['e_teff'] = 100.
    data['e_logg'] = 0.25
    data['e_fe_h'] = 0.1
    data['e_alpha_fe'] = 0.1  # NEED TO CHECK
    print 'CHECK GALAH ALPHAFE ERR'
    data['ehrv'] = 0.6
    data['fe_h'] = salaris_mh(data.fe_h, data.alpha_fe)
    data['e_fe_h'] = salaris_mh_err(
        data.fe_h, data.alpha_fe, data.e_fe_h, data.e_alpha_fe)
    cm = crossmatch_2MASS(data.ra.values, data.dec.values)
    for i in cm.keys():
        data[i] = cm[i].values
    # data = remove_duplicates(data)
    data['mag_use'] = [np.array(['J', 'H', 'K']) for i in range(len(data))]
    data['mass'] = 0.
    data['mass_error'] = -1.
    data['rho_Tg'] = 0.
    data['rho_TZ'] = 0.
    data['rho_gZ'] = 0.
    return data


def load_data():
    data = read_fits_to_pandas('/data/jls/galah/galah_dr2.fits')

    fltr = data.flag_cannon != 0.
    for i in ['teff', 'logg', 'fe_h']:
        data.loc[fltr, i] = np.nan
        data.loc[fltr, 'e_' + i] = np.nan

    fltr = data.alpha_fe == data.alpha_fe
    data.loc[fltr, 'fe_h'] = salaris_mh(data.fe_h, data.alpha_fe)
    data.loc[fltr, 'e_fe_h'] = salaris_mh_err(
        data.fe_h, data.alpha_fe, data.e_fe_h, data.e_alpha_fe)

    col_dict = {'raj2000': 'ra', 'dej2000': 'dec',
		'rv_obst': 'hrv',
                'jmag': 'J', 'hmag': 'H', 'kmag': 'K',
                'e_jmag': 'eJ', 'e_hmag': 'eH', 'e_kmag': 'eK'}
    data = data.rename(index=str, columns=col_dict)

    data['mag_use'] = [np.array(['J', 'H', 'K']) for i in range(len(data))]
    data['mass'] = 0.
    data['mass_error'] = -1.
    data['rho_Tg'] = 0.
    data['rho_TZ'] = 0.
    data['rho_gZ'] = 0.

    return data


def load_and_match(output_file='/data/jls/GaiaDR2/spectro/GALAH_input.hdf5',
                   use_dr1=False):

    loaded, data = check_and_load(output_file, 'GALAH DR2')
    if loaded:
        return data

    data = load_data()
    data = cross_match.crossmatch_gaia_spectro(data, dr1=use_dr1, epoch=2000.)

    write_input_file(data, output_file, 'GALAH DR2')
    return data
