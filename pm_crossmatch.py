from cross_match import crossmatch_gaia_spectro
import pandas as pd
import numpy as np
from spectro_data.apogee import load_data as load_apogee
from spectro_data.galah import load_data as load_galah
from spectro_data.ges import load_data as load_ges
from spectro_data.lamost import load_data as load_lamost
from spectro_data.rave import load_data as load_rave
from spectro_data.lamost import format_columns

def find_crossmatch_failures(survey):

    r = pd.read_hdf('/data/jls/GaiaDR2/spectro/%s_input.hdf5' % survey)

#    columns = ['parallax', 'parallax_error',
#               'pmra', 'pmdec',
#               'pmra_error', 'pmdec_error',
#               'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr',
#               'source_id', 'G', 'GBP', 'GRP', 'eG', 'eGRP', 'eGBP',
#               'a_g_val', 'a_g_percentile_lower', 'a_g_percentile_upper']

 #   r = r.drop(columns, axis=1)

    if survey=='GES':
        rx = load_ges()    
    if survey=='APOGEE':
        rx = load_apogee()    
    if survey=='RAVE':
        rx = load_rave(use_dr5=True)    
    if survey=='RAVEON':
        rx = load_rave(use_dr5=False)    
    if survey=='GALAH':
        rx = load_galah()    
    if survey=='LAMOST':
        rx = load_lamost()    
        rx = format_columns(rx)

    rx = crossmatch_gaia_spectro(rx, no_proper_motion=False, dist_max=5.)
    print len(rx), np.count_nonzero(r.source_id>0),np.count_nonzero(rx.source_id>0), np.count_nonzero(r.source_id!=rx.source_id)
    rx = rx[(rx.source_id > 0) & (r.source_id !=
                                   rx.source_id)].reset_index(drop=True)

    rx.to_hdf('/data/jls/GaiaDR2/spectro/%s_input_pm.hdf5' % survey, 'data')


if __name__ == '__main__':
    #find_crossmatch_failures('GES')
    #find_crossmatch_failures('RAVE')
    #find_crossmatch_failures('APOGEE')
    find_crossmatch_failures('LAMOST')
    #find_crossmatch_failures('GALAH')
    #find_crossmatch_failures('RAVEON')
