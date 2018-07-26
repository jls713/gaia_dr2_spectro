from spectro_data.apogee import load_and_match as load_apogee
from spectro_data.test_apogeetgas import load_and_match as load_apogeetgas
from spectro_data.galah import load_and_match as load_galah
from spectro_data.ges import load_and_match as load_ges
from spectro_data.lamost import load_and_match as load_lamost
from spectro_data.rave import load_and_match as load_rave
from spectro_data.segue import load_and_match as load_segue
from distances import *
from reprocess_nogaia_match import *
from spectro_data.utils import crossmatch_PS1

DR1 = False
npool = 88
random_seed = 52
with_parallax = True

name_str = '_withPRIOR_BROAD_NEW_SETUP'

with open('config.json') as config:
    output_folder = json.load(config)['dir']['output_folder']
    output_folder = '/local/scratch_2/jls/'

def flag_errors(data, covar_correction=True):
    for col in ['dist_err', 'par_err', 'dm_err',
                'log10_teff_err', 'logg_err', 'log10_av_err',
                'log10_age_err',
                'mass_err',
                'Z_err', 'dm_log10age_corr', 'log10age_Z_corr', 'dm_Z_corr']:
        data.loc[(data['flag'] == 0) & (data[col] != data[col]), 'flag'] = 6
        data.loc[(data['flag'] == 0) & (data[col] == 0.), 'flag'] = 6
    data.loc[(data['log10_age'] < -1)&(data['logg']>3.5)&(data['log10_teff']<3.9), 'flag'] = 7
    for col in ['dm_log10age_corr', 'log10age_Z_corr', 'dm_Z_corr']:
        data.loc[(data['flag'] == 0) & (np.abs(data[col]) > 1.), 'flag'] = 6
    return data

# RAVE

def run_rave_dr5():
    rave = load_rave(output_file=output_folder + 'RAVE_input.hdf5', use_dr1=DR1, use_dr5=True)
    rave = check_photometry(rave)
    run_distance_pipeline(rave,
                          output_folder + 'RAVE_DR5_distances%s.hdf5' % name_str,
                          'raveid', 'RAVE DR5',
                          npool=npool,
                          with_parallax=with_parallax,
                          thin_mag=0.1)

def run_rave_on():
    rave = load_rave(output_folder + 'RAVEON_input.hdf5', use_dr1=DR1)
    rave = check_photometry(rave)
    run_distance_pipeline(rave,
                          output_folder + 'RAVE_Cannon_distances%s.hdf5' % name_str,
                          'raveid', 'RAVE Cannon',
                          npool=npool,
                          with_parallax=with_parallax)

def run_apogee():
    apogee = load_apogee(output_file=output_folder + 'APOGEE_input.hdf5', use_dr1=DR1)
    apogee = check_photometry(apogee)
    run_distance_pipeline(apogee,
                          output_folder + 'APOGEE_distances%s.hdf5' % name_str,
                          'APOGEE_ID', 'APOGEE',
                          npool=npool,
                          with_parallax=with_parallax, thin_mag=0.1)

def run_galah():
    galah = load_galah(output_file=output_folder + 'GALAH_input.hdf5', use_dr1=DR1)
    galah = check_photometry(galah)
    run_distance_pipeline(galah,
                          output_folder + 'GALAH_distances%s.hdf5' % name_str,
                          'sobject_id', 'GALAH',
                          npool=npool,
                          with_parallax=with_parallax)

def run_ges_dr3():
    ges = load_ges(output_file=output_folder + 'GES_input.hdf5', use_dr1=DR1)
    return
    fltr = ges['mag_use'].apply(lambda x: x[0]) == 'J'
    ges.loc[fltr] = check_photometry(ges[fltr])
    run_distance_pipeline(ges,
                          output_folder + 'GES_DR3_distances%s.hdf5' % name_str,
                          'CNAME', 'GES',
                          npool=npool,
                          with_parallax=with_parallax)

def run_segue():
    segue = load_segue(output_file=output_folder+'SEGUE_input.hdf5', use_dr1=DR1) 
    segue = segue.reset_index(drop=True)
    run_distance_pipeline(segue,
                          output_folder +
                          'SEGUE_distances%s.hdf5' % (
                              name_str),
                          'specobjid', 'SEGUE',
                          npool=npool,
                          with_parallax=with_parallax)

def run_lamost():
    lamost = load_lamost(output_file=output_folder+'LAMOST_input.hdf5', use_dr1=DR1)
    lamost = check_photometry(lamost)
    STARTNUMBER, ENDNUMBER =int(sys.argv[2]), int(sys.argv[3])
    lamost = lamost.iloc[STARTNUMBER:ENDNUMBER].reset_index(drop=True)
    run_distance_pipeline(lamost,
                          output_folder +
                          'LAMOST_distances_%s_%i_%i.hdf5' % (
                              name_str, STARTNUMBER, ENDNUMBER),
                          'obsid', 'LAMOST',
                          npool=npool,
                          with_parallax=with_parallax)

def run_lamost_ps1():
    lamost = load_lamost(output_file=output_folder+'LAMOST_input.hdf5', use_dr1=DR1)
    d = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_input_PS1__withPRIOR.hdf5', 'data')
    d = d[['obsid','gP','rP','iP','egP','erP','eiP','mag_use']]
    sys_error_floor=0.01 ## From Magnier et al. (2016) and Schlafly et al. (2012)
    for i in ['egP','erP','eiP']:
        d[i]=np.sqrt(d[i]**2+sys_error_floor**2)
    lamost = lamost.drop(['mag_use'], axis=1)
    lamost = pd.merge(d, lamost, on='obsid', how='left')
    lamost = lamost.reset_index(drop=True)
    run_distance_pipeline(lamost,
                          output_folder +
                          'LAMOST_distances_PS1_%s.hdf5' % (
                              name_str),
                          'obsid', 'LAMOST',
                          npool=npool,
                          with_parallax=with_parallax)

def run_find_gridding(survey, output_name, unique_id):
    pp = pd.read_hdf(output_folder+'%s_input.hdf5'%survey)
    oo = pd.read_hdf(output_folder+'%s_distances_withPRIOR_BROAD.hdf5'%output_name)
    oo = flag_errors(oo)
    pp = pp.reset_index(drop=True)
    oo = oo.reset_index(drop=True)
    pp = pp[(oo.flag==1)|(oo.flag==6)|(oo.flag==4)].reset_index(drop=True)
    if survey!='SEGUE':
        pp = check_photometry(pp)
    run_distance_pipeline(pp,
                          output_folder + '%s_distances_FINEGRID_%s.hdf5' % (output_name, name_str),
                          unique_id, survey,
                          npool=npool,
                          with_parallax=with_parallax, mid=0.01)
    
def run_find_gridding_lamost():
    survey='LAMOST'
    output_name='LAMOST'
    unique_id='obsid'
    STARTNUMBER, ENDNUMBER =int(sys.argv[3]), int(sys.argv[4])
    pp = pd.read_hdf(output_folder+'%s_input.hdf5'%survey)
    oo = pd.read_hdf(output_folder+'%s_distances__withPRIOR_BROAD_%i_%i.hdf5'%(output_name,STARTNUMBER,ENDNUMBER))
    oo = flag_errors(oo)
    pp = pp.reset_index(drop=True)
    oo = oo.reset_index(drop=True)
    pp = pp.merge(oo[['obsid']],how='right',on='obsid').reset_index(drop=True)
    pp = pp[(oo.flag==1)|(oo.flag==6)|(oo.flag==4)].reset_index(drop=True)
    pp = check_photometry(pp)
    run_distance_pipeline(pp,
                          output_folder + '%s_distances_FINEGRID_%s_%i_%i.hdf5' % (output_name, name_str,STARTNUMBER,ENDNUMBER),
                          unique_id, survey,
                          npool=npool,
                          with_parallax=with_parallax, mid=0.01)

def run_find_gridding_ps1():
    pp = pd.read_hdf(output_folder+'LAMOST_input_PS1__withPRIOR.hdf5')
    oo = pd.read_hdf(output_folder+'LAMOST_distances_PS1__withPRIOR_BROAD.hdf5')
    oo = flag_errors(oo)
    pp = pp[(oo.flag==1)|(oo.flag==6)|(oo.flag==4)].reset_index(drop=True)
    sys_error_floor=0.01 ## From Magnier et al. (2016) and Schlafly et al. (2012)
    for i in ['egP','erP','eiP']:
        pp[i]=np.sqrt(pp[i]**2+sys_error_floor**2)
    run_distance_pipeline(pp,
                          output_folder + 'LAMOST_distances_PS1_FINEGRID_withPRIOR_BROAD.hdf5',
                          'obsid', 'LAMOST',
                          npool=npool,
                          with_parallax=with_parallax, mid=0.01)


def run_fine_gridding():
    if sys.argv[2]=='GALAH':
        run_find_gridding('GALAH', 'GALAH', 'sobject_id')
    if sys.argv[2]=='RAVEDR5':
        run_find_gridding('RAVE', 'RAVE_DR5', 'raveid')
    if sys.argv[2]=='RAVEON':
        run_find_gridding('RAVEON', 'RAVE_Cannon', 'raveid')
    if sys.argv[2]=='APOGEE':
        run_find_gridding('APOGEE', 'APOGEE', 'APOGEE_ID')
    if sys.argv[2]=='GES':
        run_find_gridding('GES', 'GES_DR3', 'CNAME')
    if sys.argv[2]=='LAMOST':
        run_find_gridding_lamost()
    if sys.argv[2]=='SEGUE':
        run_find_gridding('SEGUE', 'SEGUE', 'specobjid')
    if sys.argv[2]=='LAMOSTPS1':
        run_find_gridding_ps1()

mmp = {
        'RAVEDR5':run_rave_dr5,
        'RAVEON':run_rave_on,
        'APOGEE':run_apogee,
        'GALAH':run_galah,
        'GES':run_ges_dr3,
        'SEGUE':run_segue,
        'LAMOST':run_lamost,
        'LAMOSTPS1':run_lamost_ps1,
        'FINE':run_fine_gridding
}

if __name__=="__main__":
    mmp[sys.argv[1]]()
