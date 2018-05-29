from utils import *
sys.path.append('../')
import cross_match


def format_columns(data):
    col_dict = {'teffadopunc': 'e_teff',
                'loggadopunc': 'e_logg',
                'fehadopunc': 'e_fe_h',
                'fehadop': 'fe_h',
                'loggadop': 'logg',
                'teffadop': 'teff',
                'elodiervfinal': 'hrv',
                'elodiervfinalerr': 'e_hrv'}
    data = data.rename(index=str, columns=col_dict)
    return data


def load_data():
    data = sqlutil.get("""select s.specobjid,ra,dec,psfmag_g as g, psfmag_r as r, """
                       """psfmag_i as i,psfmagerr_g as eg, psfmagerr_r as er, """
                       """psfmagerr_i as ei, teffadop, loggadop, fehadop, teffadopunc, """
                       """loggadopunc, fehadopunc, elodiervfinal, elodiervfinalerr """
                       """from sdssdr12.sppparams as spp, """
                       """sdssdr12.specphotoall as s where spp.specobjid=s.specobjid and spp.scienceprimary=1 """
                       """and spp.programname like '%segue%' """
                       """and (spp.zwarning=0 or spp.zwarning=16) """
                       """and spp.teffadop>0. and spp.loggadop>0. and spp.fehadop>-50. and spp.flag='nnnnn'""",
                       host='cappc127', user='jason_sanders',
                       password=wsdbpassword,
                       preamb='set enable_seqscan to off; ' +
                       'set enable_mergejoin to off; ' +
                       'set enable_hashjoin to off;', asDict=True,
                       strLength=30)
    df = pd.DataFrame(columns=data.keys())
    for k in data.keys():
        df[k] = data[k]

    df = format_columns(df)

    df['mag_use'] = [np.array(['g', 'r', 'i', 'G']) for i in range(len(df))]

    df['rho_Tg'] = 0.
    df['rho_TZ'] = 0.
    df['rho_gZ'] = 0.

    df['mass'] = 0.
    df['mass_error'] = -1.

    # External errors from Lee et al. (2008a)
    df['e_teff']=np.sqrt(df['e_teff']**2+141**2)
    df['e_logg']=np.sqrt(df['e_logg']**2+0.23**2)
    df['e_fe_h']=np.sqrt(df['e_fe_h']**2+0.23**2)

    return df


def load_and_match(output_file='/data/jls/GaiaDR2/spectro/SEGUE_input.hdf5',
                   use_dr1=False):

    loaded, data = check_and_load(
        output_file, 'SEGUE 1,2,Faint catalogue')
    if loaded:
        return data

    data = load_data()
    data = cross_match.crossmatch_gaia_spectro(data, dr1=use_dr1, epoch=2000.)
    rx = cross_match.crossmatch_gaia_spectro(data, no_proper_motion=False, dist_max=5.)
    fltr = (rx.source_id > 0) & (data.source_id != rx.source_id)
    if np.count_nonzero(fltr)>0:
        data.loc[fltr]=rx.loc[fltr]

    fltr = (data.source_id<0.)
    data.loc[fltr,'mag_use']=data.loc[fltr].applymap(lambda x: np.array(['g','r','i']))

    write_input_file(data, output_file, 'SEGUE 1,2,FAINT catalogue')

    return data
