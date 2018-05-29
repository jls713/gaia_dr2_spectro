from utils import *
sys.path.append('../')
import cross_match


def grab_all_raveon():
    data = sqlutil.get("""select * from rave_on_1_0.main as r""",
                       host='cappc127', user='jason_sanders',
                       password=wsdbpassword,
                       preamb='set enable_seqscan to off; ' +
                              'set enable_mergejoin to off; ' +
                              'set enable_hashjoin to off;', asDict=True,
                       strLength=30)
    df = pd.DataFrame(columns=data.keys())
    for k in data.keys():
        df[k] = data[k]
    return df


def metallicity_from_fe_alpha(data):
    ''' inverse variance weighting '''
    D2 = 1. / data.e_o_h**2 + 1. / data.e_mg_h**2 + \
        1. / data.e_si_h**2 + 1. / data.e_ca_h**2
    yhat = data.o_h / data.e_o_h**2 + data.mg_h / data.e_mg_h**2 + \
        data.si_h / data.e_si_h**2 + data.ca_h / data.e_ca_h**2
    afe, afe_err = yhat / D2, 1. / np.sqrt(D2)
    return salaris_mh(data.fe_h, afe), \
        salaris_mh_err(data.fe_h, afe, data.e_fe_h, afe_err)


def load_data(use_dr5=False, **args):
    '''
        if use_dr5 we use the RAVE DR5 results (Teff, logg, Z) with no
        correlations between the parameters.
    '''

    data = grab_all_raveon()
    Rsort = data.sort_values(
        by='snr', ascending=False).reset_index(drop=True)
    data_radecmag = Rsort[['raveid']]
    data = Rsort[~data_radecmag.duplicated()].sort_index()
    data = data.reset_index(drop=True)

    # RAVE DR4 Teff_SPARV not necessary here
    # rave_dr4 = pd.read_csv('/data/jls/RAVEdata/RAVE_DR4.csv.gz',
    #                       na_values='\N')
    #dr4_cols = ['RAVE_OBS_ID', 'Teff_SPARV']
    # data = data.merge(rave_dr4[dr4_cols], left_on='rave_obs_id',
    #                  right_on='RAVE_OBS_ID', how='inner')

    rave_dr5 = pd.read_csv('/data/jls/RAVEdata/RAVE_DR5.csv.gz')

    dr5_cols = ['RAVE_OBS_ID',
                'Teff_N_K', 'eTeff_K',
                'logg_N_K', 'elogg_K',
                'Met_N_K', 'eMet_K',
                'Jmag_2MASS', 'eJmag_2MASS',
                'Hmag_2MASS', 'eHmag_2MASS',
                'Kmag_2MASS', 'eKmag_2MASS']

    data = data.merge(rave_dr5[dr5_cols], left_on='rave_obs_id',
                      right_on='RAVE_OBS_ID', how='inner')

    col_dict = {'Jmag_2MASS': 'J', 'eJmag_2MASS': 'eJ',
                'Hmag_2MASS': 'H', 'eHmag_2MASS': 'eH',
                'Kmag_2MASS': 'K', 'eKmag_2MASS': 'eK'}

    data = data.rename(index=str, columns=col_dict)

    # Grab correlation coefficients
    covar_pkl = '/data/gaia-eso/arc/rave-data-files/RAVE-on-covariances.pkl'
    with open(covar_pkl) as f:
        corr_data = pickle.load(f)
    df = np.vstack((corr_data[0].T,
                    np.reshape(corr_data[2][:, :3, :3],
                               (len(corr_data[2]), 9))[:, [1, 2, 5]].T)).T
    corr_cols = ['RAVE_OBS_ID', 'rho_Tg', 'rho_TZ', 'rho_gZ']
    df = pd.DataFrame(df, columns=corr_cols)
    for c in corr_cols[1:]:
        df[c] = df[c].astype(np.float64)
    data = data.merge(df, left_on='rave_obs_id',
                      right_on='RAVE_OBS_ID', how='inner')

    met, err_met = metallicity_from_fe_alpha(data)
    fltr = (met == met)
    data.loc[fltr, 'fe_h'], data.loc[fltr, 'e_fe_h'] = met[fltr], err_met[fltr]

    if use_dr5:
        data['teff'] = data['Teff_N_K']
        data['e_teff'] = data['eTeff_K']
        data['logg'] = data['logg_N_K']
        data['e_logg'] = data['elogg_K']
        data['fe_h'] = data['Met_N_K']
        data['e_fe_h'] = data['eMet_K']
        data['rho_Tg'] = 0.
        data['rho_TZ'] = 0.
        data['rho_gZ'] = 0.

    data = data.rename(index=str, columns={'ehrv': 'e_hrv'})

    data['mag_use'] = [np.array(['J', 'H', 'K', 'G'])
                       for i in range(len(data))]

    data['mass'] = 0.
    data['mass_error'] = -1.

    return data


def load_and_match(output_file='/data/jls/GaiaDR2/spectro/RAVE_input.hdf5',
                   use_dr5=False, use_dr1=False):

    name = 'Cannon'
    if use_dr5:
        name = 'DR5'

    loaded, data = check_and_load(output_file, 'RAVE %s' % name)
    if loaded:
        return data

    data = load_data(use_dr5=use_dr5)
    print 'Loaded. Now xmatch.'
    data = cross_match.crossmatch_gaia_spectro(data, dr1=use_dr1, epoch=2000.)

    write_input_file(data, output_file, 'RAVE %s' % name)

    return data
