import pandas as pd
from astropy.table import Table, vstack, join
import numpy as np
import actions

print 'NEED A COVARIANCE FLAG -- SOME STARS COVAR IS CORRECT'

def flag_errors(data, covar_correction=True):
    for col in ['dist_err', 'par_err', 'dm_err',
                'log10_teff_err', 'logg_err', 'log10_av_err',
                'log10_age_err',
                'mass_err',
                'Z_err', 'dm_log10age_corr', 'log10age_Z_corr', 'dm_Z_corr']:
        data['flag'][(data['flag'] == 0) & (data[col] != data[col])] = 6
        data['flag'][(data['flag'] == 0) & (data[col] == 0.)] = 6
        #data[col][(data['flag']==0)&(data[col]==0.)]=\
	#	np.nanpercentile(data[col][data['flag']==0],1.)
    data['flag'][(data['log10_age'] < -1)&(data['logg']>3.5)&(data['log10_teff']<3.9)] = 7
    #for col in ['R_err']:
    #    data['flag'][(data['flag'] == 0) & (data[col] != data[col])] = 7

    #### data['source_id'][data['source_id'] < 0.] = np.nan

    # Fix correlation error in distances.py
    if covar_correction:
        data['log10age_Z_corr']*=data['log10_age']*data['Z']/data['log10_age_err']/data['Z_err']
        data['dm_Z_corr']*=data['dm']*data['Z']/data['dm_err']/data['Z_err']
        data['dm_log10age_corr']*=data['dm']*data['log10_age']/data['log10_age_err']/data['dm_err']
    for col in ['dm_log10age_corr', 'log10age_Z_corr', 'dm_Z_corr']:
        data['flag'][(data['flag'] == 0) & (np.abs(data[col]) > 1.)] = 6
    return data


extra_meta = [
              'flag=6 nan errors from pipeline (likely only one isochrone point available)',
              'flag=7 pre-main-sequence star with age < 100 Myr (likely binary star?)',
              'duplicated -- 1 if duplicate Gaia source_id -- we keep sources in the order: '
	      'APOGEE, GALAH, GES, RAVE-ON, RAVE and LAMOST',
              'best -- 1 if flag=0, duplicated=0 and valid Gaia source_id'
             ]


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


def build_master_output():

    # LAMOST
    la = Table.read('/data/jls/GaiaDR2/spectro/LAMOST_distances.hdf5')
    la = la[:1589000]
    la2 = Table.read('/data/jls/GaiaDR2/spectro/LAMOST_distances_2.hdf5')
    lm = vstack([la, la2])
    # lfill = Table.read('/data/jls/GaiaDR2/spectro/LAMOST_distances_fillin.hdf5')
    # lm = fillin(lm, lfill, 'obsid')
    lm = flag_errors(lm)
    lm.meta['COMMENT'] = la.meta['COMMENT'] + extra_meta
    lm.write('/data/jls/GaiaDR2/spectro/LAMOST_distances_MASTER.hdf5',
             path='data', format='hdf5',
             compression=True,
             serialize_meta=True,
             overwrite=True)

    # RAVE
    r = Table.read('/data/jls/GaiaDR2/spectro/RAVE_DR5_distances.hdf5')
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/RAVE_DR5_distances_fillin.hdf5')
    r = fillin(r, rfill, 'raveid')

    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/RAVE_DR5_distances_extra.hdf5')
    r = fillin(r, rfill, 'raveid')

    r = flag_errors(r)
    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/RAVE_DR5_distances_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)

    # RAVE-ON
    r = Table.read('/data/jls/GaiaDR2/spectro/RAVE_Cannon_distances.hdf5')
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/RAVE_Cannon_distances_fillin.hdf5')
    r = fillin(r, rfill, 'raveid')

    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/RAVE_Cannon_distances_extra.hdf5')
    r = fillin(r, rfill, 'raveid')

    r = flag_errors(r)
    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/RAVE_Cannon_distances_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)

    # APOGEE
    r = Table.read('/data/jls/GaiaDR2/spectro/APOGEE_distances.hdf5')

    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/APOGEE_distances_fillin.hdf5')
    r = fillin(r, rfill, 'APOGEE_ID')

    # rfill = Table.read(
    #     '/data/jls/GaiaDR2/spectro/APOGEE_distances_extra.hdf5')
    # r = fillin(r, rfill, 'APOGEE_ID')

    r = flag_errors(r)
    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/APOGEE_distances_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)

    # GALAH
    r = Table.read('/data/jls/GaiaDR2/spectro/GALAH_distances.hdf5')
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/GALAH_distances_fillin.hdf5')
    r = fillin(r, rfill, 'sobject_id')

    # rfill = Table.read(
    #     '/data/jls/GaiaDR2/spectro/GALAH_distances_extra.hdf5')
    # r = fillin(r, rfill, 'sobject_id')

    r = flag_errors(r)
    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/GALAH_distances_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)

    # GES
    r = Table.read('/data/jls/GaiaDR2/spectro/GES_distances.hdf5')
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/GES_distances_fillin.hdf5')
    r = fillin(r, rfill, 'CNAME')

    # rfill = Table.read(
    #     '/data/jls/GaiaDR2/spectro/GES_distances_extra.hdf5')
    # r = fillin(r, rfill, 'CNAME')

    r = flag_errors(r)
    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/GES_distances_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)

def test_segue():
    actions.init_()
    # SEGUE
    r = Table.read('/data/jls/GaiaDR2/spectro/SEGUE_distances_withPRIOR.hdf5')
    r = flag_errors(r, covar_correction=False)
    #rfill = Table.read(
    #    '/data/jls/GaiaDR2/spectro/SEGUE_distances_FINEGRID__withPRIOR.hdf5')
    #rfill = flag_errors(rfill, covar_correction=False)
    #r = fillin(r, rfill, 'specobjid')

    linput = pd.read_hdf('/data/jls/GaiaDR2/spectro/SEGUE_input_MASTER.hdf5')
    r=actions.process_actions(linput[:4],r[:4],None)
    r.write('/data/jls/GaiaDR2/spectro/tmp.hdf5',
             path='data', format='hdf5',
             compression=True,
             serialize_meta=True,
             overwrite=True)


def build_master_output_withPRIOR():
    actions.init_()

    # LAMOST
    fls = ['0_500000', '500000_1000000', '900000_1000000','1000000_1500000', '1350000_1500000',
           '1500000_2000000', '1700000_2000000','2000000_2500000', '2100000_2500000',
           '2500000_3500000', '2567800_3500000']
    lm = Table()
    for f in fls:
        la = Table.read('/data/jls/GaiaDR2/spectro/LAMOST_distances__withPRIOR_%s.hdf5'%f)
        la = la[la['l']==la['l']]
        lm = vstack([lm,la])
    lfill = Table.read('/data/jls/GaiaDR2/spectro/LAMOST_distances_pm__withPRIOR.hdf5')
    lm = fillin(lm, lfill, 'obsid')
    lm = flag_errors(lm)
    lfill = Table.read('/data/jls/GaiaDR2/spectro/LAMOST_distances_PS1__withPRIOR.hdf5')
    lfill = flag_errors(lfill, covar_correction=False)
    lm = fillin(lm, lfill, 'obsid')
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/LAMOST_distances_FINEGRID__withPRIOR.hdf5')
    rfill = flag_errors(rfill, covar_correction=False)
    lm = fillin(lm, rfill, 'obsid')
    linput = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_input_MASTER.hdf5')
    lm=actions.process_actions(linput,lm,None)
    lm.meta['COMMENT'] = la.meta['COMMENT'] + extra_meta
    lm.write('/data/jls/GaiaDR2/spectro/LAMOST_distances_withPRIOR_MASTER.hdf5',
             path='data', format='hdf5',
             compression=True,
             serialize_meta=True,
             overwrite=True)
    lm['survey']='LAMOST'
    del linput

    # RAVE
    r = Table.read('/data/jls/GaiaDR2/spectro/RAVE_DR5_distances_withPRIOR.hdf5')
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/RAVE_DR5_distances_pm_withPRIOR.hdf5')
    r = fillin(r, rfill, 'raveid')
    r = flag_errors(r)
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/RAVE_DR5_distances_FINEGRID__withPRIOR.hdf5')
    rfill = flag_errors(rfill, covar_correction=False)
    r = fillin(r, rfill, 'raveid')
    linput = pd.read_hdf('/data/jls/GaiaDR2/spectro/RAVE_input_MASTER.hdf5')
    r=actions.process_actions(linput,r,None)
    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/RAVE_DR5_distances_withPRIOR_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)
    r['survey']='RAVEDR5'
    meta = r.meta['COMMENT'][3]
    r.meta['COMMENT']=None
    lm = vstack([lm,r])
    lm.meta['COMMENT'].insert(3,meta)
    del r
    del linput


    # RAVE-ON
    r = Table.read('/data/jls/GaiaDR2/spectro/RAVE_Cannon_distances_withPRIOR.hdf5')
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/RAVE_Cannon_distances_pm_withPRIOR.hdf5')
    r = fillin(r, rfill, 'raveid')
    r = flag_errors(r)
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/RAVE_Cannon_distances_FINEGRID__withPRIOR.hdf5')
    rfill = flag_errors(rfill, covar_correction=False)
    r = fillin(r, rfill, 'raveid')

    linput = pd.read_hdf('/data/jls/GaiaDR2/spectro/RAVEON_input_MASTER.hdf5')
    r=actions.process_actions(linput,r,None)
    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/RAVE_Cannon_distances_withPRIOR_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)
    r['survey']='RAVEON'
    meta = r.meta['COMMENT'][3]
    r.meta['COMMENT']=None
    lm = vstack([lm,r])
    lm.meta['COMMENT'].insert(3,meta)
    del r
    del linput


    # APOGEE
    r = Table.read('/data/jls/GaiaDR2/spectro/APOGEE_distances_withPRIOR.hdf5')

    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/APOGEE_distances_pm_withPRIOR.hdf5')
    r = fillin(r, rfill, 'APOGEE_ID')
    r = flag_errors(r)

    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/APOGEE_distances_FINEGRID__withPRIOR.hdf5')
    rfill = flag_errors(rfill, covar_correction=False)
    r = fillin(r, rfill, 'APOGEE_ID')

    linput = pd.read_hdf('/data/jls/GaiaDR2/spectro/APOGEE_input_MASTER.hdf5')
    r=actions.process_actions(linput,r,None)
    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/APOGEE_distances_withPRIOR_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)
    r['survey']='APOGEE'
    meta = r.meta['COMMENT'][3]
    r.meta['COMMENT']=None
    lm = vstack([lm,r])
    lm.meta['COMMENT'].insert(3,meta)
    del r
    del linput


    # GALAH
    r = Table.read('/data/jls/GaiaDR2/spectro/GALAH_distances_withPRIOR.hdf5')
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/GALAH_distances_pm_withPRIOR.hdf5')
    r = fillin(r, rfill, 'sobject_id')
    r = flag_errors(r)
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/GALAH_distances_GIANTS_withPRIOR.hdf5')
    rfill = flag_errors(rfill, covar_correction=False)
    r = fillin(r, rfill, 'sobject_id')

    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/GALAH_distances_FINEGRID__withPRIOR.hdf5')
    rfill = flag_errors(rfill, covar_correction=False)
    r = fillin(r, rfill, 'sobject_id')

    linput = pd.read_hdf('/data/jls/GaiaDR2/spectro/GALAH_input_MASTER.hdf5')
    r=actions.process_actions(linput,r,None)
    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/GALAH_distances_withPRIOR_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)
    r['survey']='GALAH'
    meta = r.meta['COMMENT'][3]
    r.meta['COMMENT']=None
    lm = vstack([lm,r])
    lm.meta['COMMENT'].insert(3,meta)
    del r
    del linput


    # GES
    r = Table.read('/data/jls/GaiaDR2/spectro/GES_DR3_distances_withPRIOR.hdf5')
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/GES_DR3_distances_pm_withPRIOR.hdf5')
    r = fillin(r, rfill, 'CNAME')
    r = flag_errors(r)
    rfill = Table.read(
        '/data/jls/GaiaDR2/spectro/GES_DR3_distances_FINEGRID__withPRIOR.hdf5')
    rfill = flag_errors(rfill, covar_correction=False)
    r = fillin(r, rfill, 'CNAME')
    linput = pd.read_hdf('/data/jls/GaiaDR2/spectro/GES_input_MASTER.hdf5')
    r=actions.process_actions(linput,r,None)
    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/GES_DR3_distances_withPRIOR_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)
    r['survey']='GES'
    meta = r.meta['COMMENT'][3]
    r.meta['COMMENT']=None
    lm = vstack([lm,r])
    lm.meta['COMMENT'].insert(3,meta)
    del r
    del linput

    # SEGUE
    r = Table.read('/data/jls/GaiaDR2/spectro/SEGUE_distances_withPRIOR.hdf5')
    r = flag_errors(r, covar_correction=False)
    rfill = Table.read(
       '/data/jls/GaiaDR2/spectro/SEGUE_distances_FINEGRID__withPRIOR.hdf5')
    rfill = flag_errors(rfill, covar_correction=False)
    r = fillin(r, rfill, 'specobjid')

    linput = pd.read_hdf('/data/jls/GaiaDR2/spectro/SEGUE_input_MASTER.hdf5')
    r=actions.process_actions(linput,r,None)

    r.meta['COMMENT'] += extra_meta
    r.write('/data/jls/GaiaDR2/spectro/SEGUE_distances_withPRIOR_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)
    r['survey']='SEGUE'
    meta = r.meta['COMMENT'][3]
    r.meta['COMMENT']=None
    lm = vstack([lm,r])
    lm.meta['COMMENT'].insert(3,meta)
    del r
    del linput

    # Keep non-duplicate entries
    lm['duplicated']=0
    lp = lm[['source_id','survey','flag','vz_err']].to_pandas()
    lp['indx']=np.arange(len(lp))
    lp.survey = pd.Categorical(lp.survey,
                      categories=["APOGEE","GALAH","GES","RAVEON","RAVE","LAMOST", "SEGUE"],
                      ordered=True)
    lp = lp.sort_values(['flag','survey','vz_err'])
    fltr = lp['source_id'].duplicated()
    lm['duplicated'][lp['indx'].values[fltr]]=1

    lm['best']=0
    lm['best'][(lm['source_id']==lm['source_id'])&(lm['flag']==0)&(lm['duplicated']==0)]=1

    lm.meta['COMMENT'][0]='Distances, ages, masses and extinctions for spectroscopic surveys combined with Gaia DR2'

    lm.write('/data/jls/GaiaDR2/spectro/distances_withPRIOR_MASTER.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)

def build_input_master():
    for s, fld in zip(['LAMOST','GES','APOGEE','RAVE','RAVEON','GALAH'],
                      ['obsid','CNAME','APOGEE_ID','raveid','raveid','sobject_id']):
        p = pd.read_hdf('/data/jls/GaiaDR2/spectro/%s_input.hdf5'%s)
        p2 = pd.read_hdf('/data/jls/GaiaDR2/spectro/%s_input_pm.hdf5'%s)
        p = p.merge(p2, on=fld, how='left', suffixes=['','_copy'])
        for col in p.columns:
            if col[-5:]=='_copy':
                p.loc[p[col]==p[col],col[:-5]]=p[col][p[col]==p[col]]
                p=p.drop([col], axis=1)
        p.to_hdf('/data/jls/GaiaDR2/spectro/%s_input_MASTER.hdf5'%s, 'data')

if __name__ == '__main__':
    #build_input_master()
    #build_master_output()
    build_master_output_withPRIOR()
    #test_segue()
