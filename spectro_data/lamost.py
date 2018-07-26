from utils import *
sys.path.append('../')
import cross_match


def format_columns(data):
    col_dict = {'teff_err': 'e_teff',
                'logg_err': 'e_logg',
                'feh_err': 'e_fe_h',
                'feh': 'fe_h',
                'rv': 'hrv',
                'rv_err': 'e_hrv'}
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

    df['mag_use'] = [np.array(['J', 'H', 'K', 'G']) for i in range(len(df))]

    df['rho_Tg'] = 0.
    df['rho_TZ'] = 0.
    df['rho_gZ'] = 0.

    output_file = '/data/jls/GaiaDR2/spectro/lamost_cannon/LAMOST_results.hdf5'
    t = pd.read_hdf(output_file)
    fltr = (t.TEFF>4000.)&(t.TEFF<5250.)&(t.LOGG>1.)&(t.LOGG<3.3)&(t.M_H>-1.5)&(t.M_H<0.5)
    fltr &= (t.r_chi_sq < 3.) #& (t.in_convex_hull == True)
    t = t[fltr].reset_index(drop=True)

    df = df.merge(t[['obsid','mass','mass_error']], how='left')
    fltr = df.mass!=df.mass

    df.loc[fltr,'mass'] = 0.
    df.loc[fltr,'mass_error'] = -1.

    return df


def load_and_match(output_file='/data/jls/GaiaDR2/spectro/LAMOST_input.hdf5',
                   use_dr1=False):

    loaded, data = check_and_load(
        output_file, 'LAMOST DR3 A, F, G, K catalogue')
    if loaded:
        return data

    data = load_data()
    data = format_columns(data)
    rx = cross_match.crossmatch_gaia_spectro(data, no_proper_motion=False, dist_max=5.)
    data = cross_match.crossmatch_gaia_spectro(data, dr1=use_dr1, epoch=2000.)
    fltr = (rx.source_id > 0) & (data.source_id != rx.source_id)
    if np.count_nonzero(fltr)>0:
        data.loc[fltr]=rx.loc[fltr]
    data = data.reset_index(drop=True)

    write_input_file(data, output_file, 'LAMOST DR3 A,F,G,K catalogue')

    return data
