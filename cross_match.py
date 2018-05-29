import numpy as np
import sys
import sqlutil
from login import wsdbpassword


def crossmatch_gaia(data, dr1=False, epoch=2000,
                    dist_max=5., no_proper_motion=True,
                    phot_g_cut=None):
    """
        Cross match list of ra,dec of stars with Gaia catalogues.
        (ra,dec) vectors of RA and Dec in deg
        epoch epoch of observation -- if epoch='ref_epoch' use m.ref_epoch
        dist_max maximum angular separation to consider (in arcsec)
        max_epoch_diff = 20 years -- each star can have separate epoch.
    """
    ra, dec = data.ra.values, data.dec.values
    gaia_catalogue = 'gaia_dr2.gaia_source'
    if isinstance(epoch, basestring):
        join_strng = 'q3c_join(m.ra,m.dec,s.ra,s.dec,%0.8f) and q3c_join_pm(s.ra,s.dec,s.pmra,s.pmdec,s.ref_epoch,' % (dist_max / 3600. * 10.) +\
            'm.ra,m.dec,%s,20.,%0.8f)' % (epoch, dist_max / 3600.)
        dist_strng = 'q3c_dist_pm(s.ra,s.dec,s.pmra,s.pmdec,s.ref_epoch,' +\
            'm.ra,m.dec,%s)' % (epoch)
    else:
        join_strng = 'q3c_join(m.ra,m.dec,s.ra,s.dec,%0.8f) and q3c_join_pm(s.ra,s.dec,s.pmra,s.pmdec,s.ref_epoch,' % (dist_max / 3600. * 10.) +\
            'm.ra,m.dec,%0.2f,20.,%0.8f)' % (epoch, dist_max / 3600.)
        dist_strng = 'q3c_dist_pm(s.ra,s.dec,s.pmra,s.pmdec,s.ref_epoch,' +\
            'm.ra,m.dec,%0.2f)' % (epoch)
    if dr1:
        gaia_catalogue = 'gaia_dr1.gaia_source'
    if no_proper_motion:
        join_strng = 'q3c_join(m.ra,m.dec,s.ra,s.dec,%0.8f)' % (
            dist_max / 3600.)
        dist_strng = 'q3c_dist(m.ra,m.dec,s.ra,s.dec)'
    phot_g_strng = ''
    if phot_g_cut is not None:
        phot_g_strng = 's.phot_g_mean_mag<%0.5f and' % phot_g_cut
    rqst = """
        select tt.* from mytable as m left join lateral
        (select *, %s as dist from %s as s where %s %s order by %s  asc limit 1)
        as tt on  true  order by xid """ % \
        (dist_strng, gaia_catalogue, phot_g_strng, join_strng, dist_strng)
    print rqst
    add_data = sqlutil.local_join(rqst,
                                  'mytable',
                                  (ra, dec, np.arange(len(dec))),
                                  ('ra', 'dec', 'xid'),
                                  host='cappc127',
                                  user='jason_sanders',
                                  password=wsdbpassword,
                                  preamb='set enable_seqscan to off; ' +
                                  'set enable_mergejoin to off; ' +
                                  'set enable_hashjoin to off;',
                                  asDict=True,
                                  strLength=30)

    columns = ['parallax', 'parallax_error',
               'pmra', 'pmdec',
               'pmra_error', 'pmdec_error',
               'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr',
               'source_id',
               'phot_g_mean_mag', 'phot_rp_mean_mag',
               'phot_bp_mean_mag', 'phot_g_mean_flux',
               'phot_rp_mean_flux', 'phot_bp_mean_flux',
               'phot_g_mean_flux_error',
               'phot_rp_mean_flux_error', 'phot_bp_mean_flux_error',
               'a_g_val', 'a_g_percentile_lower', 'a_g_percentile_upper',
               'dist']

    for c in add_data.keys():
        if c in columns:
            data[c] = add_data[c]

    phot_keys = ['g', 'rp', 'bp']

    if dr1:
        phot_keys = ['g']

    for c in phot_keys:
        flx_err = data['phot_%s_mean_flux_error' % c]
        flx_err /= data['phot_%s_mean_flux' % c]
        data['ephot_%s_mean_mag' % c] = 2.5 / np.log(10.) * flx_err
        data = data.drop(['phot_%s_mean_flux' % c,
                          'phot_%s_mean_flux_error' % c], axis=1)

    col_dict = {'phot_g_mean_mag': 'G',
                'phot_rp_mean_mag': 'GRP',
                'phot_bp_mean_mag': 'GBP',
                'ephot_g_mean_mag': 'eG',
                'ephot_rp_mean_mag': 'eGRP',
                'ephot_bp_mean_mag': 'eGBP'}

    data = data.rename(index=str, columns=col_dict)

    return data


def crossmatch_gaia_ids(data):
    """
        Cross match list of ra,dec of stars with Gaia catalogues.
        (ra,dec) vectors of RA and Dec in deg
        epoch epoch of observation -- if epoch='ref_epoch' use m.ref_epoch
        dist_max maximum angular separation to consider (in arcsec)
        max_epoch_diff = 20 years -- each star can have separate epoch.
    """
    gaia_catalogue = 'gaia_dr2.gaia_source'
    rqst = """
        select tt.* from %s as tt, mytable as m where tt.source_id=m.source_id""" % \
        (gaia_catalogue)
    print rqst
    add_data = sqlutil.local_join(rqst,
                                  'mytable',
                                  (data.source_id.values,),
                                  ('source_id',),
                                  host='cappc127',
                                  user='jason_sanders',
                                  password=wsdbpassword,
                                  preamb='set enable_seqscan to off; ' +
                                  'set enable_mergejoin to off; ' +
                                  'set enable_hashjoin to off;',
                                  asDict=True,
                                  strLength=30)

    columns = ['parallax', 'parallax_error',
               'pmra', 'pmdec',
               'pmra_error', 'pmdec_error',
               'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr',
               'source_id',
               'phot_g_mean_mag', 'phot_rp_mean_mag',
               'phot_bp_mean_mag', 'phot_g_mean_flux',
               'phot_rp_mean_flux', 'phot_bp_mean_flux',
               'phot_g_mean_flux_error',
               'phot_rp_mean_flux_error', 'phot_bp_mean_flux_error',
               'a_g_val', 'a_g_percentile_lower', 'a_g_percentile_upper']

    for c in add_data.keys():
        if c in columns:
            data[c] = add_data[c]

    phot_keys = ['g', 'rp', 'bp']

    for c in phot_keys:
        flx_err = data['phot_%s_mean_flux_error' % c]
        flx_err /= data['phot_%s_mean_flux' % c]
        data['ephot_%s_mean_mag' % c] = 2.5 / np.log(10.) * flx_err
        data = data.drop(['phot_%s_mean_flux' % c,
                          'phot_%s_mean_flux_error' % c], axis=1)

    col_dict = {'phot_g_mean_mag': 'G',
                'phot_rp_mean_mag': 'GRP',
                'phot_bp_mean_mag': 'GBP',
                'ephot_g_mean_mag': 'eG',
                'ephot_rp_mean_mag': 'eGRP',
                'ephot_bp_mean_mag': 'eGBP'}

    data = data.rename(index=str, columns=col_dict)

    return data


def crossmatch_gaia_spectro_id(data):
    df = crossmatch_gaia_ids(data)

    sys_error_floor = 0.02  # Approximate systematic errors of 20mmag in G
    # Also reflects the degree to which the isochrones are good.
    phot = ['G', 'GBP', 'GRP']
    for c in phot:
        df['e%s' % c] = np.sqrt(
            df['e%s' % c]**2 + sys_error_floor**2)
    fltr = (df.parallax != df.parallax)
    df.loc[fltr, 'parallax'] = 0.
    df.loc[fltr, 'parallax_error'] = -1.
    return df


def crossmatch_gaia_spectro(data, dr1=False, epoch=2000, dist_max=5.,
                            no_proper_motion=True, phot_g_cut=None):
    df = crossmatch_gaia(data, dr1=dr1, epoch=epoch,
                         dist_max=dist_max,
                         no_proper_motion=no_proper_motion,
                         phot_g_cut=phot_g_cut)

    sys_error_floor = 0.02  # Approximate systematic errors of 20mmag in G
    # Also reflects the degree to which the isochrones are good.
    phot = ['G', 'GBP', 'GRP']
    if dr1:
        phot = ['G']
    for c in phot:
        df['e%s' % c] = np.sqrt(
            df['e%s' % c]**2 + sys_error_floor**2)
    fltr = (df.parallax != df.parallax)
    df.loc[fltr, 'parallax'] = 0.
    df.loc[fltr, 'parallax_error'] = -1.
    return df


import pandas as pd


def crossmatch_2mass_gaia_sourceid(data):

    gaia_catalogue = 'gaia_dr2_aux.gaia_source_2mass_xm'
    rqst = """
        select * from mytable as m left join %s as s on s.source_id=m.source_id""" % \
        (gaia_catalogue)
    print rqst
    add_data = sqlutil.local_join(rqst,
                                  'mytable',
                                  (data.source_id.values, ),
                                  ('source_id',),
                                  host='cappc127',
                                  user='jason_sanders',
                                  password=wsdbpassword,
                                  preamb='set enable_seqscan to off; ' +
                                  'set enable_mergejoin to off; ' +
                                  'set enable_hashjoin to off;',
                                  asDict=True,
                                  strLength=30)
    return pd.DataFrame.from_dict(add_data)


def crossmatch_2mass_ids_gaia(data, dist_max=2.):

    ra, dec = data.ra.values, data.dec.values
    gaia_catalogue = 'gaia_dr2_aux.gaia_source_2mass_xm'
    join_strng = 'q3c_join(m.ra,m.dec,s.t_ra,s.t_decl,%0.8f)' % (
        dist_max / 3600.)
    dist_strng = 'q3c_dist(m.ra,m.dec,s.t_ra,s.t_decl)'
    rqst = """
        select tt.* from mytable as m left join lateral
        (select * from %s as s where %s order by %s  asc limit 1)
        as tt on  true  order by xid """ % \
        (gaia_catalogue, join_strng, dist_strng)
    print rqst
    add_data = sqlutil.local_join(rqst,
                                  'mytable',
                                  (ra, dec, np.arange(len(dec))),
                                  ('ra', 'dec', 'xid'),
                                  host='cappc127',
                                  user='jason_sanders',
                                  password=wsdbpassword,
                                  preamb='set enable_seqscan to off; ' +
                                  'set enable_mergejoin to off; ' +
                                  'set enable_hashjoin to off;',
                                  asDict=True,
                                  strLength=30)
    return pd.DataFrame.from_dict(add_data)


def crossmatch_2mass_ids_gaia_version2(data, dist_max=10.):

    ra, dec = data.ra.values, data.dec.values
    J = data.J.values
    join_strng = 'q3c_join(m.ra,m.dec,s.ra,s.dec,%0.8f)' % (
        dist_max / 3600.)
    dist_strng = 'q3c_dist(m.ra,m.dec,s.ra,s.dec)'
    rqst = """
        select tt.* from mytable as m left join lateral
        (select * from gaia_dr2_aux.gaia_source_2mass_xm as s where %s
        and abs(s.j_m-m.j)<1e-3
        order by %s  asc limit 1)
        as tt on true order by xid""" % (join_strng, dist_strng)
    # and abs(s.j_m-m.j)<1e-3
    print rqst
    add_data = sqlutil.local_join(rqst,
                                  'mytable',
                                  (ra, dec, J, np.arange(len(dec))),
                                  ('ra', 'dec', 'j', 'xid'),
                                  host='cappc127',
                                  user='jason_sanders',
                                  password=wsdbpassword,
                                  preamb='set enable_seqscan to on; ' +
                                  'set enable_mergejoin to on; ' +
                                  'set enable_hashjoin to on;',
                                  asDict=True,
                                  strLength=30)
    return pd.DataFrame.from_dict(add_data)
