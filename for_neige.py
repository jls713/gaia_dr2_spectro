# Run with flat age prior

from spectro_data.apogee import load_and_match as load_apogee
from distances import *
from reprocess_nogaia_match import *
from spectro_data.utils import crossmatch_PS1
from astropy.table import Table, join

DR1 = False
npool = 88
random_seed = 52
with_parallax = True

name_str = '_withPRIOR_FLAT_081018'

with open('config.json') as config:
    output_folder = json.load(config)['dir']['output_folder']
    output_folder = '/local/scratch_2/jls/'

def flag_errors(data):
    for col in ['dist_err', 'par_err', 'dm_err',
                'log10_teff_err', 'logg_err', 'log10_av_err',
                'log10_age_err',
                'mass_err',
                'Z_err', 'dm_log10age_corr', 'log10age_Z_corr', 'dm_Z_corr']:
        data['flag'][(data['flag'] == 0) & (data[col] != data[col])] = 6
        data['flag'][(data['flag'] == 0) & (data[col] == 0.)] = 6
    data['flag'][(data['log10_age'] < -1)&(data['logg']>3.5)&(data['log10_teff']<3.9)] = 7
    for col in ['dm_log10age_corr', 'log10age_Z_corr', 'dm_Z_corr']:
        data['flag'][(data['flag'] == 0) & (np.abs(data[col]) > 1.)] = 6
    return data

def run_apogee():
    apogee = load_apogee(output_file=output_folder + 'APOGEE_input.hdf5', use_dr1=DR1)
    apogee = check_photometry(apogee)
    run_distance_pipeline(apogee,
                          output_folder + 'APOGEE_distances%s.hdf5' % name_str,
                          'APOGEE_ID', 'APOGEE',
                          npool=npool,
                          with_parallax=with_parallax, thin_mag=0.1)

def run_find_gridding(survey, output_name, unique_id):
    pp = pd.read_hdf(output_folder+'%s_input.hdf5'%survey)
    oo = pd.read_hdf(output_folder+'%s_distances%s.hdf5'%(output_name,name_str))
    oo = flag_errors(oo)
    pp = pp.reset_index(drop=True)
    oo = oo.reset_index(drop=True)
    pp = pp[(oo.flag==1)|(oo.flag==6)|(oo.flag==4)].reset_index(drop=True)
    run_distance_pipeline(pp,
                          output_folder + '%s_distances_FINEGRID_%s.hdf5' % (output_name, name_str),
                          unique_id, survey,
                          npool=npool,
                          with_parallax=with_parallax, mid=0.01)

def fillin(r, rfill, fld):
    r['indx'] = np.arange(len(r))
    rfill['indx'] = np.arange(len(rfill))
    r2 = join(r[[fld, 'indx']],
              rfill[[fld, 'indx']],
              keys=fld, join_type='left')
    r[r2['indx_1'][~r2['indx_2'].mask]] = \
        rfill[r2['indx_2'][~r2['indx_2'].mask]]
    r.remove_column('indx')
    return r

def combine_files(output_name='APOGEE'):
    pp = Table.read(output_folder+'%s_distances%s.hdf5'%(output_name,name_str))
    pp = flag_errors(pp)
    pp2 = Table.read(output_folder+'%s_distances_FINEGRID_%s.hdf5'%(output_name,name_str))
    pp2 = flag_errors(pp2)
    pp = fillin(pp,pp2,'APOGEE_ID')
    pp.meta['COMMENT'] = pp.meta['COMMENT'] + ['Parallax zeropoint -0.029mas included',
					       'Flat age prior',
					       'systematic G error of 0.01 mag added'] 
    pp.write('/data/jls/GaiaDR2/spectro/APOGEE_distances_flatageprior_parzero_sysG_forneige.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)

if __name__=="__main__":
    #run_apogee()
    #run_find_gridding('APOGEE', 'APOGEE', 'APOGEE_ID')
    combine_files()
