import pandas as pd
from astropy.table import Table, vstack, join
from astropy import units as u
import numpy as np
import actions
from reprocess_nogaia_match import check_photometry

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
    lfill = Table.read('/data/jls/GaiaDR2/spectro/LAMOST_distances_PS1_FINEGRID_withPRIOR.hdf5')
    lfill = flag_errors(lfill, covar_correction=False)
    lm = fillin(lm, lfill, 'obsid')
    lfill = Table.read('/data/jls/GaiaDR2/spectro/LAMOST_distances_PS1_FINEGRID_SATURATED_withPRIOR.hdf5')
    lfill = flag_errors(lfill, covar_correction=False)
    lm = fillin(lm, lfill, 'obsid')

    linput = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_input_MASTER.hdf5')

    ### Adjust LAMOST radial velocities -- add 4.5
    linput['hrv']+=4.5
    ##########

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
    rfill = Table.read(
       '/data/jls/GaiaDR2/spectro/SEGUE_distances_withPRIOR_parallax_fillin.hdf5')
    rfill = flag_errors(rfill, covar_correction=False)
    r = fillin(r, rfill, 'specobjid')
    #############################
    ## Failed on SEGUE, no parallax -- adding here
    rfill = Table.read('/data/jls/GaiaDR2/spectro/SEGUE_distances_withPRIOR_parallax_fillin.hdf5')
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
        p = Table.from_pandas(pd.read_hdf('/data/jls/GaiaDR2/spectro/%s_input.hdf5'%s))
        p2 = Table.from_pandas(pd.read_hdf('/data/jls/GaiaDR2/spectro/%s_input_pm.hdf5'%s))
        columns = [fld, 'parallax', 'parallax_error',
               'pmra', 'pmdec',
               'pmra_error', 'pmdec_error',
               'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr',
               'source_id',
               'G','eG','GBP','eGBP','GRP','eGRP',
               'a_g_val', 'a_g_percentile_lower', 'a_g_percentile_upper',
               'dist']
        p['indx']=np.arange(len(p))
        p = join(p, p2[columns], keys=fld, join_type='left', uniq_col_name='{col_name}{table_name}',
		 table_names=['','_copy'])
        fltr = ~p['source_id_copy'].mask
        p['source_id'][fltr]=p['source_id_copy'][fltr]
        del p['source_id_copy']
        for col in p.colnames:
            if col[-5:]=='_copy':
                p[col[:-5]][fltr]=p[col][fltr]
                del p[col]
        p['mag_use']=[','.join(pp) for pp in p['mag_use']]
        for i in p.colnames:
            if p[i].dtype==type(object):
                p[i]=p[i].astype('|S')
        p.sort('indx')
        del p['indx']
        p.write('/data/jls/GaiaDR2/spectro/%s_input_MASTER.hdf5'%s, path='data',
                format='hdf5',compression=True,serialize_meta=True,overwrite=True)


def fix_master_file():
    d = Table.read('/data/jls/GaiaDR2/spectro/distances_withPRIOR_MASTER.hdf5')
    #############################
    # Keep non-duplicate entries
    d['duplicated']=0
    lp = d[['source_id','survey','flag','vz_err']].to_pandas()
    lp['indx']=np.arange(len(lp))
    lp.survey = pd.Categorical(lp.survey,
                      categories=["APOGEE","GALAH","GES","RAVEON","RAVE","LAMOST", "SEGUE"],
                      ordered=True)
    lp = lp.sort_values(['flag','survey','vz_err'])
    fltr = lp['source_id'].duplicated()
    d['duplicated'][lp['indx'].values[fltr]]=1

    d['best']=0
    d['best'][(d['source_id']==d['source_id'])&(d['flag']==0)&(d['duplicated']==0)]=1
    ##############################


    del d['source_id']
    d['source_id']=-1
    d['ra']=0.*u.deg
    d['dec']=0.*u.deg
    d['angular_separation']=np.nan*u.arcsec
    d['mag_use']=''
    d['mag_use']=d['mag_use'].astype('|S10')

    for s, fld in zip(['LAMOST','GES','APOGEE','RAVE','RAVEON','GALAH','SEGUE'],
                      ['obsid','CNAME','APOGEE_ID','raveid','raveid','sobject_id','specobjid']):

        p = pd.read_hdf('/data/jls/GaiaDR2/spectro/%s_input.hdf5'%s)

        if s=='LAMOST':
            p = check_photometry(p)
            pp = pd.read_hdf('/data/jls/GaiaDR2/spectro/LAMOST_input_PS1__withPRIOR.hdf5')
            saturated = (pp.gP<13.5)|(pp.rP<13.5)|(pp.iP<13.5)
            pp.loc[saturated,'mag_use']=pp.applymap(lambda x: np.array(['G','GBP','GRP']))
            p = pd.merge(p, pp[['obsid','mag_use']], on='obsid', how='left')
            p['mag_use']=p['mag_use_x']
            p.loc[p['mag_use_y']==p['mag_use_y'],'mag_use']=p['mag_use_y'][p['mag_use_y']==p['mag_use_y']]
            del p['mag_use_y']
            del p['mag_use_x']
        elif s=='GES':
            fltr = p['mag_use'].apply(lambda x: x[0]) == 'J'
            p.loc[fltr] = check_photometry(p[fltr])
        elif s=='SEGUE':
            pass
        else:
            p = check_photometry(p)

        p['mag_use']=[','.join(pp) for pp in p['mag_use']]
        p['mag_use']=p['mag_use'].astype(str)

        p = Table.from_pandas(p)

        if s!='SEGUE':
            p2 = Table.from_pandas(pd.read_hdf('/data/jls/GaiaDR2/spectro/%s_input_pm.hdf5'%s))
            columns = [fld,
               'source_id',
               'ra', 'dec']
            p = join(p, p2[columns], keys=fld, join_type='left', uniq_col_name='{col_name}{table_name}', table_names=['','_copy'])
            #p['dist'][(p['source_id_copy']>=0)&(~p['source_id_copy'].mask)]=p['dist_copy'][(p['source_id_copy']>=0)&(~p['source_id_copy'].mask)]
            p['source_id'][(p['source_id_copy']>=0)&(~p['source_id_copy'].mask)]=p['source_id_copy'][(p['source_id_copy']>=0)&(~p['source_id_copy'].mask)]
            del p['source_id_copy']
            #del p['dist_copy']
            del p['ra_copy']
            del p['dec_copy']
        p_i = pd.read_hdf('/local/scratch_2/jls/%s_input.hdf5'%s)
        p['angular_separation']=p_i['dist'].values*3600.*u.arcsec
        #p['dist'].name='dist_crossmatch'
        if s!='RAVE':
            survey=s
        else:
            survey='RAVEDR5'
        p['survey']=survey
        print len(p), np.count_nonzero(p['ra'])
        p = p[['source_id',fld,'survey','ra','dec','angular_separation','mag_use']]
        d = join(d,p,join_type='left',keys=[fld,'survey'], uniq_col_name='{col_name}{table_name}',table_names=['','_2'])
        for f in ['ra','dec','angular_separation','mag_use','source_id']:
            d[f][(~d['source_id_2'].mask)]=d[f+'_2'][(~d['source_id_2'].mask)]
            del d[f+'_2']
        print np.count_nonzero(d['ra'])

    d['sobject_id']=d['sobject_id'].astype(np.int64)
    d['specobjid']=d['specobjid'].astype(np.int64)

    lm = d[['survey','flag','vz_err']].to_pandas()
    lm['source_id']=d['source_id']
    # Keep non-duplicate entries
    lm['duplicated']=0
    lp = lm.copy()
    lp['indx']=np.arange(len(lp))
    lp.survey = pd.Categorical(lp.survey,
                      categories=["APOGEE","GALAH","GES","RAVEON","RAVE","LAMOST", "SEGUE"],
                      ordered=True)
    lp = lp.sort_values(['flag','survey','vz_err'])
    fltr = lp['source_id'].duplicated()
    lm.loc[lp['indx'].values[fltr],'duplicated']=1

    lm['best']=0
    lm.loc[(lm['source_id']>=0)&(lm['flag']==0)&(lm['duplicated']==0),'best']=1

    d['duplicated']=lm['duplicated']
    d['best']=lm['best']
    d.meta['COMMENT'][3]='survey is a string giving the name of the survey the star was taken from -- APOGEE, GALAH, GES, RAVEON, RAVEDR5, LAMOST, SEGUE'

    d.meta['COMMENT'][10]='ra,dec are the equatorial coordinates from the spectroscopic surveys (in degrees) -- we use the spectroscopic survey ra,dec instead of Gaia ra,dec so cross-match failures with Gaia still have values.'
    d.meta['COMMENT'][11]='angular_separation is the cross-match distance in arcsec between the spectroscopic catalogue and Gaia DR2 (if the proper motion has been considered, it is the epoch correct distance).'
    d.meta['COMMENT'][12]='mag_use is a comma-separated string giving the magnitudes used in the pipeline. J,H,K are 2MASS bands, G,GBP,GRP are Gaia bands, gP,rP,iP are Pan-STARRS bands, g,r,i are SDSS bands, Jv,Hv,Kv are VISTA bands'

    d.write('/data/jls/GaiaDR2/spectro/distances_withPRIOR_fixedint64.hdf5',
            path='data', format='hdf5',
            compression=True,
            serialize_meta=True,
            overwrite=True)


if __name__ == '__main__':
    #build_input_master()
    #build_master_output()
    #build_master_output_withPRIOR()
    #test_segue()
    fix_master_file()
