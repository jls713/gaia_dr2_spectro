import numpy as np
import sys
import sqlutil
from login import wsdbpassword

def crossmatch_gaia(data, dr1=False, epoch=2000, dist_max=5.):
    """
        Cross match list of ra,dec of stars with Gaia catalogues.
        (ra,dec) vectors of RA and Dec in deg
        epoch epoch of observation
        dist_max maximum angular separation to consider (in arcsec)
        max_epoch_diff = 15 years
        <-- Not totally sure what this means, ask Sergey!
    """
    ra, dec = data.ra.values, data.dec.values
    gaia_catalogue = 'gaia_dr2.gaia_source'
    join_strng = 'q3c_join_pm(s.ra,s.dec,s.pmra,s.pmdec,s.ref_epoch,' +\
        'm.ra,m.dec,%0.2f,15.,%0.8f)' % (epoch, dist_max / 3600.)
    dist_strng = 'q3c_dist_pm(s.ra,s.dec,s.pmra,s.pmdec,s.ref_epoch,' +\
        'm.ra,m.dec,%0.2f)' % (epoch)
    if dr1:
        gaia_catalogue = 'gaia_dr1.gaia_source'
        join_strng = 'q3c_join(m.ra,m.dec,s.ra,s.dec,%0.8f)' % (
            dist_max / 3600.)
        dist_strng = 'q3c_dist(m.ra,m.dec,s.ra,s.dec)'
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

    columns = ['parallax', 'parallax_error',
               'pmra', 'pmdec',
               'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr',
               'source_id',
               'phot_g_mean_mag', 'phot_rp_mean_mag',
               'phot_bp_mean_mag', 'phot_g_mean_flux',
               'phot_rp_mean_flux', 'phot_bp_mean_flux',
               'phot_g_mean_flux_error',
               'phot_rp_mean_flux_error', 'phot_bp_mean_flux_error']

    for c in add_data.keys():
        if c in columns:
            data[c] = add_data[c]

    phot_keys = ['g', 'rp', 'bp']

    if dr1:
        phot_keys = ['g']

    sys_error_floor = 0.02  # Approximate systematic errors of 20mmag in G
    # Also reflects the degree to which the isochrones are good.

    for c in phot_keys:
        flx_err = data['phot_%s_mean_flux_error' % c]
        flx_err /= data['phot_%s_mean_flux' % c]
        data['ephot_%s_mean_mag' % c] = 2.5 / np.log(10.) * flx_err
        data['ephot_%s_mean_mag' % c] = np.sqrt(
            data['ephot_%s_mean_mag' % c]**2 + sys_error_floor**2)
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


def crossmatch_gaia_spectro(data, dr1=False, epoch=2000, dist_max=5.):
    df = crossmatch_gaia(data, dr1=dr1, epoch=epoch, dist_max=dist_max)
    fltr = (df.parallax != df.parallax)
    df.loc[fltr, 'parallax'] = 0.
    df.loc[fltr, 'parallax_error'] = -1.
    return df
