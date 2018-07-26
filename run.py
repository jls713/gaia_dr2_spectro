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
TEST = False #False
TESTNUMBER = 100  # 100
STARTNUMBER = 1500000
ENDNUMBER = 2000000
npool = 88
random_seed = 52
with_parallax = True
name_str = ''
if not with_parallax:
    name_str = '_preG'
if TEST:
    name_str = '_TEST'

name_str = '_withPRIOR'

with open('config.json') as config:
    output_folder = json.load(config)['dir']['output_folder']
    output_folder = '/local/scratch_1/jls/'
# RAVE


def run_rave_on():
    rave = load_rave(output_folder + 'RAVEON_input.hdf5', use_dr1=DR1)
    if TEST:
        rave = rave.sample(n=TESTNUMBER,
                           random_state=random_seed).reset_index(drop=True)
    rave = check_photometry(rave)
    run_distance_pipeline(rave,
                          output_folder + 'RAVE_Cannon_distances%s.hdf5' % name_str,
                          'raveid', 'RAVE Cannon',
                          npool=npool,
                          with_parallax=with_parallax)


def run_rave_on_fillin():
    rave = load_rave(output_folder + 'RAVEON_input.hdf5', use_dr1=DR1)
    rave2 = Table.read(output_folder + 'RAVE_Cannon_distances.hdf5')
    rave = check_photometry(rave[rave2['flag'] == 3].reset_index(drop=True))
    run_distance_pipeline(rave,
                          output_folder + 'RAVE_Cannon_distances_fillin.hdf5',
                          'raveid', 'RAVE Cannon',
                          npool=npool,
                          with_parallax=with_parallax)


def run_rave_on_extra():
    rave = load_rave(output_folder + 'RAVEON_input_pm.hdf5', use_dr1=DR1)
    if TEST:
        rave = rave.sample(n=TESTNUMBER,
                           random_state=random_seed).reset_index(drop=True)
    rave = check_photometry(rave)
    run_distance_pipeline(rave,
                          output_folder + 'RAVE_Cannon_distances_pm%s.hdf5' % name_str,
                          'raveid', 'RAVE Cannon',
                          npool=npool,
                          with_parallax=with_parallax)


def run_rave_dr5():
    rave = load_rave(use_dr1=DR1, use_dr5=True)
    if TEST:
        rave = rave.sample(n=TESTNUMBER,
                           random_state=random_seed).reset_index(drop=True)
    rave = check_photometry(rave)
    run_distance_pipeline(rave,
                          output_folder + 'RAVE_DR5_distances%s.hdf5' % name_str,
                          'raveid', 'RAVE DR5',
                          npool=npool,
                          with_parallax=with_parallax)


def run_rave_dr5_fillin():
    rave = load_rave(use_dr1=DR1, use_dr5=True)
    rave2 = Table.read(output_folder + 'RAVE_DR5_distances.hdf5')
    rave = check_photometry(rave[rave2['flag'] == 3].reset_index(drop=True))
    run_distance_pipeline(rave,
                          output_folder + 'RAVE_DR5_distances_fillin.hdf5',
                          'raveid', 'RAVE DR5',
                          npool=npool,
                          with_parallax=with_parallax)


def run_rave_dr5_extra():
    rave = load_rave(output_folder + 'RAVE_input_pm.hdf5',
                     use_dr1=DR1, use_dr5=True)
    if TEST:
        rave = rave.sample(n=TESTNUMBER,
                           random_state=random_seed).reset_index(drop=True)
    rave = check_photometry(rave)
    run_distance_pipeline(rave,
                          output_folder + 'RAVE_DR5_distances_pm%s.hdf5' % name_str,
                          'raveid', 'RAVE DR5',
                          npool=npool,
                          with_parallax=with_parallax)
# Gaia-ESO


def run_ges():
    ges = load_ges(use_dr1=DR1)
    if TEST:
        ges = ges.sample(n=TESTNUMBER,
                         random_state=random_seed).reset_index(drop=True)
    fltr = ges['mag_use'].apply(lambda x: x[0]) == 'J'
    ges.loc[fltr] = check_photometry(ges[fltr])
    run_distance_pipeline(ges,
                          output_folder + 'GES_distances%s.hdf5' % name_str,
                          'CNAME', 'GES',
                          npool=npool,
                          with_parallax=with_parallax)


def run_ges_dr3():
    ges = load_ges(use_dr1=DR1)
    if TEST:
        ges = ges.sample(n=TESTNUMBER,
                         random_state=random_seed).reset_index(drop=True)
    fltr = ges['mag_use'].apply(lambda x: x[0]) == 'J'
    ges.loc[fltr] = check_photometry(ges[fltr])
    run_distance_pipeline(ges,
                          output_folder + 'GES_DR3_distances%s.hdf5' % name_str,
                          'CNAME', 'GES',
                          npool=npool,
                          with_parallax=with_parallax)


def run_ges_dr3_extra():
    ges = pd.read_hdf('/data/jls/GaiaDR2/spectro/GES_input_pm.hdf5')
    fltr = ges['mag_use'].apply(lambda x: x[0]) == 'J'
    ges.loc[fltr] = check_photometry(ges[fltr])
    run_distance_pipeline(ges,
                          output_folder + 'GES_DR3_distances_pm%s.hdf5' % name_str,
                          'CNAME', 'GES',
                          npool=npool,
                          with_parallax=with_parallax)


def run_ges_fillin():
    ges = load_ges(use_dr1=DR1)
    ges2 = Table.read(output_folder + 'GES_distances.hdf5')
    ges = check_photometry(ges[ges2['flag'] == 3].reset_index(drop=True))
    run_distance_pipeline(ges,
                          output_folder + 'GES_distances_fillin.hdf5',
                          'CNAME', 'GES',
                          npool=npool,
                          with_parallax=with_parallax)

# APOGEE
# APOGEE


def run_apogee():
    apogee = load_apogee(use_dr1=DR1)
    if TEST:
        apogee = apogee.sample(n=TESTNUMBER,
                               random_state=random_seed).reset_index(drop=True)
    apogee = check_photometry(apogee)
    run_distance_pipeline(apogee,
                          output_folder + 'APOGEE_distances%s.hdf5' % name_str,
                          'APOGEE_ID', 'APOGEE',
                          npool=npool,
                          with_parallax=with_parallax)


def run_apogee_extra():
    apogee = pd.read_hdf('/data/jls/GaiaDR2/spectro/APOGEE_input_pm.hdf5')
    apogee = check_photometry(apogee)
    run_distance_pipeline(apogee,
                          output_folder + 'APOGEE_distances_pm%s.hdf5' % name_str,
                          'APOGEE_ID', 'APOGEE',
                          npool=npool,
                          with_parallax=with_parallax)


def run_apogee_fillin():
    apogee = load_apogee(use_dr1=DR1)
    apogee2 = Table.read(output_folder + 'APOGEE_distances.hdf5')
    apogee = check_photometry(
        apogee[apogee2['flag'] == 3].reset_index(drop=True))
    run_distance_pipeline(apogee,
                          output_folder + 'APOGEE_distances_fillin.hdf5',
                          'APOGEE_ID', 'APOGEE',
                          npool=npool,
                          with_parallax=with_parallax)


# SEGUE

def run_segue():
    segue = load_segue(use_dr1=DR1) 
    if TEST:
        segue = segue.sample(n=TESTNUMBER,
                               random_state=random_seed).reset_index(drop=True)
    segue = segue.reset_index(drop=True)
    run_distance_pipeline(segue,
                          output_folder +
                          'SEGUE_distances_TEST_%s.hdf5' % (
                              name_str),
                          'specobjid', 'SEGUE',
                          npool=npool,
                          with_parallax=with_parallax)

# LAMOST


def run_lamost():
    lamost = load_lamost(use_dr1=DR1)
    if TEST:
        lamost = lamost.sample(n=TESTNUMBER,
                               random_state=random_seed).reset_index(drop=True)
    lamost = check_photometry(lamost)
    lamost = lamost.iloc[STARTNUMBER:ENDNUMBER].reset_index(drop=True)
    run_distance_pipeline(lamost,
                          output_folder +
                          'LAMOST_distances_%s_%i_%i.hdf5' % (
                              name_str, STARTNUMBER, ENDNUMBER),
                          'obsid', 'LAMOST',
                          npool=npool,
                          with_parallax=with_parallax)


def run_lamost_extra():
    lamost = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_input_pm.hdf5')
    col_dict = {'teff_err': 'e_teff',
                'logg_err': 'e_logg',
                'feh_err': 'e_fe_h',
                'feh': 'fe_h',
                'rv': 'hrv',
                'rv_err': 'e_hrv'}
    lamost = lamost.rename(index=str, columns=col_dict)
    lamost = check_photometry(lamost)
    run_distance_pipeline(lamost,
                          output_folder +
                          'LAMOST_distances_pm_%s.hdf5' % (name_str),
                          'obsid', 'LAMOST',
                          npool=npool,
                          with_parallax=with_parallax)
# GALAH


def run_galah():
    galah = load_galah(use_dr1=DR1)
    if TEST:
        galah = galah.sample(n=TESTNUMBER,
                             random_state=random_seed).reset_index(drop=True)
    galah = check_photometry(galah)
    run_distance_pipeline(galah,
                          output_folder + 'GALAH_distances%s.hdf5' % name_str,
                          'sobject_id', 'GALAH',
                          npool=npool,
                          with_parallax=with_parallax)


def run_galah_extra():
    galah = pd.read_hdf('/data/jls/GaiaDR2/spectro/GALAH_input_pm.hdf5')
    galah = check_photometry(galah)
    run_distance_pipeline(galah,
                          output_folder + 'GALAH_distances_pm%s.hdf5' % name_str,
                          'sobject_id', 'GALAH',
                          npool=npool,
                          with_parallax=with_parallax)


def run_galah_fillin():
    galah = load_galah(use_dr1=DR1)
    galah2 = Table.read(output_folder + 'GALAH_distances.hdf5')
    galah = check_photometry(galah[galah2['flag'] == 3].reset_index(drop=True))
    run_distance_pipeline(galah,
                          output_folder + 'GALAH_distances_fillin.hdf5',
                          'sobject_id', 'GALAH',
                          npool=npool,
                          with_parallax=with_parallax)

def run_galah_giants_mass():
    galah = pd.read_hdf('/data/jls/GaiaDR2/spectro/GALAH_input_MASTER.hdf5')
    payel = pd.read_hdf('/data/jls/GaiaDR2/spectro/GALAH_input_payel_masses.hdf5')
    payel = payel[payel.massflag==0].reset_index(drop=True)
    galah = galah.drop(['mass','mass_error'], axis=1)
    galah = pd.merge(galah, payel[['sobject_id','mass','mass_error']], 
                     on='sobject_id', how='inner')
    galah = check_photometry(galah)
    run_distance_pipeline(galah,
                          output_folder + 'GALAH_distances_GIANTS%s.hdf5' % name_str,
                          'sobject_id', 'GALAH',
                          npool=npool,
                          with_parallax=with_parallax)
     

# APOGEE TGAS TEST


def run_apogeetgas_test():
    apogee = load_apogeetgas()
    print len(apogee)
    if TEST:
        apogee = apogee.sample(n=TESTNUMBER,
                               random_state=random_seed).reset_index(drop=True)
    print len(apogee)
    run_distance_pipeline(apogee,
                          output_folder +
                          'APOGEETGAS_TEST_distances_random100.hdf5',
                          'location_id', 'APOGEE TGAS',
                          npool=npool,
                          with_parallax=with_parallax)


def run_apogeetgas_test_withGaia():
    apogee = load_apogeetgas()
    apogee['eG'] = 0.025
    apogee['e_hrv'] = apogee['evlos']
    apogee['mag_use'] = \
        apogee.applymap(
            lambda x: np.array(['J', 'H', 'K', 'G']))['mag_use']
    if TEST:
        apogee = apogee.sample(n=TESTNUMBER,
                               random_state=random_seed).reset_index(drop=True)
    run_distance_pipeline(apogee,
                          output_folder +
                          'APOGEETGAS_TEST_distances_random100_G.hdf5',
                          'location_id', 'APOGEE TGAS',
                          npool=npool)

def run_lamost_ps1():
    dd = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_input_MASTER.hdf5')
    rdd = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_distances_withPRIOR_MASTER.hdf5')
    dd = dd.loc[rdd.flag.values==3].reset_index(drop=True)
    ps1 = crossmatch_PS1(dd.ra.values,dd.dec.values)
    for k in ps1.keys():
        if k!='ra' and k!='dec':
            dd[k] = ps1[k]
    
    col_dict = {'gpsfmagerr': 'egP',
                'rpsfmagerr': 'erP',
                'ipsfmagerr': 'eiP',
                'gpsfmag': 'gP',
                'rpsfmag': 'rP',
                'ipsfmag': 'iP'}
    dd = dd.rename(index=str, columns=col_dict)
    dd['mag_use'] = \
        dd.applymap(
            lambda x: np.array(['gP', 'rP', 'iP', 'G']))['mag_use']
    dd = dd.reset_index(drop=True)
    fltr = (dd.gP!=dd.gP)
    dd.loc[fltr,'mag_use']=dd.loc[fltr].applymap(lambda x: np.array(['G','GBP','GRP']))
    dd.to_hdf(output_folder + 'LAMOST_input_PS1_%s.hdf5' % name_str, 'data')
    if TEST:
        dd = dd.sample(n=TESTNUMBER,
                               random_state=random_seed).reset_index(drop=True)
    run_distance_pipeline(dd,
                          output_folder + 'LAMOST_distances_PS1_%s.hdf5' % (name_str),
                          'obsid', 'LAMOST',
                          npool=npool,
                          with_parallax=with_parallax, mid=None)
 


def run_find_gridding(survey, output_name, unique_id):
    pp = pd.read_hdf('/data/jls/GaiaDR2/spectro/%s_input_MASTER.hdf5'%survey)
    oo = pd.read_hdf('/data/jls/GaiaDR2/spectro/%s_distances_withPRIOR_MASTER.hdf5'%output_name)
    pp = pp.reset_index(drop=True)
    oo = oo.reset_index(drop=True)
    pp = pp[(oo.flag==1)|(oo.flag==6)].reset_index(drop=True)
    if survey=='GALAH':
        payel = pd.read_hdf('/data/jls/GaiaDR2/spectro/GALAH_input_payel_masses.hdf5')
        payel = payel[payel.massflag==0].reset_index(drop=True)
        pp = pp.drop(['mass','mass_error'], axis=1)
        pp = pd.merge(pp, payel[['sobject_id','mass','mass_error']], 
                     on='sobject_id', how='left')
        fltr = (pp['mass']!=pp['mass'])
        pp.loc[fltr,'mass']=0.
        pp.loc[fltr,'mass_error']=0.
    if survey!='SEGUE':
        pp = check_photometry(pp)
    run_distance_pipeline(pp,
                          output_folder + '%s_distances_FINEGRID_%s.hdf5' % (output_name, name_str),
                          unique_id, survey,
                          npool=npool,
                          with_parallax=with_parallax, mid=0.01)

def run_parallax_fix_segue():
    pp = pd.read_hdf('/data/jls/GaiaDR2/spectro/SEGUE_input_MASTER.hdf5')    
    oo = pd.read_hdf('/data/jls/GaiaDR2/spectro/SEGUE_distances_withPRIOR_MASTER.hdf5')    
    pp = pp.reset_index(drop=True)
    oo = oo.reset_index(drop=True)
    pp = pp[(oo.flag==4)].reset_index(drop=True)
    pp['parallax']=0.
    pp['parallax_error']=-1.
    run_distance_pipeline(pp,
                          output_folder + 'SEGUE_distances_withPRIOR_parallax_fillin.hdf5',
                          'specobjid', 'SEGUE',
                          npool=npool,
                          with_parallax=with_parallax, mid=None, thin_mag=-1.)
                          
def run_find_gridding_ps1():
    pp = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_input_PS1__withPRIOR.hdf5')
    oo = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_distances_PS1__withPRIOR.hdf5')
    pp = pp[(oo.flag==1)|(oo.flag==6)].reset_index(drop=True)
    run_distance_pipeline(pp,
                          output_folder + 'LAMOST_distances_PS1_FINEGRID_withPRIOR.hdf5',
                          'obsid', 'LAMOST',
                          npool=npool,
                          with_parallax=with_parallax, mid=0.01)

if __name__=="__main__":
    run_parallax_fix_segue()
    # run_segue()
    #run_lamost_ps1()
    #run_find_gridding('GALAH', 'GALAH', 'sobject_id')
    #run_find_gridding('RAVE', 'RAVE_DR5', 'raveid')
    #run_find_gridding('RAVEON', 'RAVE_Cannon', 'raveid')
    #run_find_gridding('APOGEE', 'APOGEE', 'APOGEE_ID')
    #run_find_gridding('GES', 'GES_DR3', 'CNAME')
    #run_find_gridding('LAMOST', 'LAMOST', 'obsid')
    #run_find_gridding('SEGUE', 'SEGUE', 'specobjid')
    #run_find_gridding_ps1()
    # run_galah_giants_mass()
    # run_rave_on_extra()
    # run_rave_on()
    # run_rave_dr5()
    # run_ges()
    # run_ges_dr3()
    # run_apogee()
    # run_lamost()
    # run_galah()
    # run_rave_on_fillin()
    # run_rave_dr5_fillin()
    # run_ges_fillin()
    # run_apogee_fillin()
    # run_lamost()
    # run_galah_fillin()
    # run_apogeetgas_test()
    # run_apogeetgas_test_withGaia()
