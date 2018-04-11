from utils import *
sys.path.append('../')
import cross_match


def format_columns(data):
    col_dict = {'teff_err': 'e_teff',
                'logg_err': 'e_logg',
                'feh_err': 'e_fe_h',
                'feh': 'fe_h',
                'rv': 'hrv'}
    data = data.rename(index=str, columns=col_dict)
    return data


def load_data():

    data = sqlutil.get("""select * from lamost_dr3.stellar as l""",
                       host='cappc127', user='jason_sanders',
                       password=wsdbpassword,
                       preamb='set enable_seqscan to off; ' +
                       'set enable_mergejoin to off; ' +
                       'set enable_hashjoin to off;', asDict=True,
                       strLength=30)
    df = pd.DataFrame(columns=data.keys())
    for k in data.keys():
        df[k] = data[k]

    cm2MASS = crossmatch_2MASS(df.ra.values, df.dec.values)
    cm2MASS = quality_2MASS_phot(cm2MASS)
    cmSDSS = crossmatch_SDSS(df.ra.values, df.dec.values)

    for k in ['J', 'eJ', 'H', 'eH', 'K', 'eK']:
        df[k] = cm2MASS[k].values
    for k in cmSDSS.keys():
        df[k] = cmSDSS[k]

    df['mag_use'] = [np.array(['J', 'H', 'K']) for i in range(len(df))]

    df['mass'] = 0.
    df['mass_error'] = -1.

    df['rho_Tg'] = 0.
    df['rho_TZ'] = 0.
    df['rho_gZ'] = 0.
    
    return df


def load_and_match(output_file='/data/jls/GaiaDR2/spectro/LAMOST_input.hdf5',
                   use_dr1=False):

    loaded, data = check_and_load(output_file, 'LAMOST DR3 A, F, G, K catalogue')
    if loaded:
        return data

    data = load_data()
    data = format_columns(data)
    data = cross_match.crossmatch_gaia_spectro(data, dr1=use_dr1, epoch=2000.)

    write_input_file(data, output_file, 'LAMOST DR3 A,F,G,K catalogue')

    return data
