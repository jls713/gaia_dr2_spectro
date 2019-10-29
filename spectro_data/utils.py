import pandas as pd
import fitsio
import numpy as np
from astropy.table import Table
from login import wsdbpassword
import os
#import cPickle as pickle
import _pickle as pickle
import sys
sys.path.append('/data/jls/astrolibpy/utils')
import sqlutil

kwargs = {'host': 'cappc127',
          'user': 'jason_sanders',
          'password': wsdbpassword, 
          'preamb': 'set enable_seqscan to off; ' +
          'set enable_mergejoin to off; ' +
          'set enable_hashjoin to off;',
          'asDict': True,
          'strLength': 30}


def salaris_mh(fe_h, a_fe):
    return fe_h + np.log10(0.638 * pow(10., a_fe) + 0.362)


def salaris_mh_err(fe_h, a_fe, fe_h_err, a_fe_err):
    da = (2.30259 * pow(10., a_fe)) / (0.567398 + pow(10., a_fe))
    return np.sqrt(fe_h_err**2 + da**2 * a_fe_err**2)


def read_fits_to_pandas(data_file):
    ''' Reads a fits file into a pandas DataFrame. The column names in the
        dataframe are set by those in the fits file '''
    output_data = pd.DataFrame()
    input_data = fitsio.FITS(data_file)[1].read().byteswap().newbyteorder()
    for i in input_data.dtype.descr:
        if(isinstance(input_data[i[0]][0], np.ndarray)):
            for c, j in enumerate(input_data[i[0]].T):
                if(isinstance(j, np.ndarray)):
                    continue
                else:
                    output_data[i[0] + str(c)] = j
        else:
            output_data[i[0]] = input_data[i[0]]
    return output_data


def check_and_load(output_file, survey):
    if os.path.isfile(output_file):
        return True,pd.read_hdf(output_file)
        #tbl = Table.read(output_file)
        tbl = Table.from_pandas(pd.read_hdf(output_file))
        tbl.meta['COMMENT'] = 'Input file for distance pipeline using the ' + \
                              'results from %s spectroscopic survey.' % survey
        tbl = tbl.to_pandas()
        return True, tbl
    else:
        return False, None


def write_input_file(data, output_file, survey):
    data.to_hdf(output_file, 'data', comp_level=9)
    tbl = Table.from_pandas(data)
    tbl.meta['COMMENT'] = 'Input file for distance pipeline using the ' + \
                          'results from %s spectroscopic survey.' % survey
    return
    tbl.write(output_file, path='data', format='hdf5',
              compression=True,
              serialize_meta=True,
              overwrite=True)


def crossmatch_2MASS(ra, dec, dist_max=5.):
    """
        Cross match list of ra,dec of stars with 2MASS.
    """
    data = sqlutil.local_join("""
        select j_m as J,j_msigcom as eJ,h_m as H,h_msigcom as eH,
               k_m as K,k_msigcom as eK, j_psfchi, h_psfchi, k_psfchi,
               ph_qual, rd_flg, cc_flg from mytable as m
        left join lateral (select * from twomass.psc as s
        where q3c_join(m.ra, m.dec,s.ra,s.decl,%0.5f/3600)
        order by q3c_dist(m.ra,m.dec,s.ra,s.decl) asc limit 1)
        as tt on  true  order by xid """ % dist_max,
                              'mytable',
                              (ra, dec, np.arange(len(dec))),
                              ('ra', 'dec', ' xid'), **kwargs)
    df2MASS = pd.DataFrame(columns=data.keys())
    for c in data.keys():
        df2MASS[c] = data[c]
    df2MASS = df2MASS.rename(index=str, columns={'j': 'J', 'h': 'H', 'k': 'K',
                                                 'ej': 'eJ', 'eh': 'eH', 'ek': 'eK'})
    return df2MASS


def crossmatch_SDSS(ra, dec, dist_max=5.):
    """
        Cross match list of ra,dec of stars with SDSS.
    """
    data = sqlutil.local_join("""
        select psfMag_g as g,psfMagErr_g as eg,
        psfMag_r as r,psfMagErr_r as er,
        psfMag_i as i,psfMagErr_i as ei from mytable as m
        left join lateral (select * from sdssdr12.specphotoall as s
        where s.type=6 and mode=1 and
        q3c_join(m.ra, m.dec,s.ra,s.dec,%0.5f/3600)
        order by q3c_dist(m.ra,m.dec,s.ra,s.dec) asc limit 1)
        as tt on  true  order by xid """ % dist_max,
                              'mytable', (ra, dec, np.arange(
                                  len(dec))), ('ra', 'dec', 'xid'),
                              host='cappc127', user='jason_sanders',
                              password=wsdbpassword,
                              preamb='set enable_seqscan to off; ' +
                              'set enable_mergejoin to off; ' +
                              'set enable_hashjoin to off;',
                              asDict=True, strLength=30)
    return data

def crossmatch_PS1(ra, dec, dist_max=5.):
    """
        Cross match list of ra,dec of stars with Pan-STARRS.
    """
    data = sqlutil.local_join("""
        select * from mytable as m
        left join lateral (select * from panstarrs_dr1.stackobjectthin as s where
        q3c_join(m.ra, m.dec,s.ra,s.dec,%0.5f/3600)
        and s.ipsfmag - s.ikronmag < 0.05
        and (s.ginfoflag3&panstarrs_dr1.detectionflags3('STACK_PRIMARY'))>0
        order by q3c_dist(m.ra,m.dec,s.ra,s.dec) asc limit 1)
        as tt on  true  order by xid """ % dist_max,
                              'mytable', (ra, dec, np.arange(
                                  len(dec))), ('ra', 'dec', 'xid'),
                              host='cappc127', user='jason_sanders',
                              password=wsdbpassword,
                              preamb='set enable_seqscan to off; ' +
                              'set enable_mergejoin to off; ' +
                              'set enable_hashjoin to off;',
                              asDict=True, strLength=30)
    return data

def quality_2MASS_phot(tmass_data):

    tmass_data['rd_flg'] = tmass_data['rd_flg'].astype(str).str.zfill(3)
    fltrJ, fltrH, fltrK = 0, 0, 0
    for a in ['A', 'B', 'C', 'D']:
        fltrJ |= tmass_data['ph_qual'].str.contains('%s..' % a, regex=True)
        fltrH |= tmass_data['ph_qual'].str.contains('.%s.' % a, regex=True)
        fltrK |= tmass_data['ph_qual'].str.contains('..%s' % a, regex=True)
    fltrJ &= (tmass_data['cc_flg'].str.contains('0..', regex=True))
    fltrH &= (tmass_data['cc_flg'].str.contains('.0.', regex=True))
    fltrK &= (tmass_data['cc_flg'].str.contains('..0', regex=True))
    fltrJ &= ~((tmass_data['j_psfchi'] > 3.) & (
        tmass_data['rd_flg'].str.contains('2..', regex=True)))
    fltrH &= ~((tmass_data['h_psfchi'] > 3.) & (
        tmass_data['rd_flg'].str.contains('.2.', regex=True)))
    fltrK &= ~((tmass_data['k_psfchi'] > 3.) & (
        tmass_data['rd_flg'].str.contains('..2', regex=True)))
    tmass_data.loc[~fltrJ, 'J'] = np.nan
    tmass_data.loc[~fltrJ, 'eJ'] = np.nan
    tmass_data.loc[~fltrH, 'H'] = np.nan
    tmass_data.loc[~fltrH, 'eH'] = np.nan
    tmass_data.loc[~fltrK, 'K'] = np.nan
    tmass_data.loc[~fltrK, 'eK'] = np.nan

    return tmass_data


def crossmatch_vista(ra, dec, dist_max=5.):
    """
        Cross match list of ra,dec of stars with VISTA.
    """
    data1 = sqlutil.local_join("""
        select japercormag4,japercormag4_err,
        hapercormag4,hapercormag4_err,
        kapercormag4,kapercormag4_err,ebv from mytable as m
        left join lateral (select * from vhs_1708.des as s
        where prim=1 and mergedclass=-1
        and jerrorbit=0 and kerrorbit=0
        and q3c_join(m.ra, m.dec,s.ra,s.dec,%0.5f/3600)
        order by q3c_dist(m.ra,m.dec,s.ra,s.dec) asc limit 1)
        as tt on  true  order by xid """ % dist_max,
                               'mytable', (ra, dec, np.arange(
                                   len(dec))), ('ra', 'dec', 'xid'), **kwargs)
    data2 = sqlutil.local_join("""
        select japercormag4,japercormag4_err,
        hapercormag4,hapercormag4_err,
        kapercormag4,kapercormag4_err,ebv from mytable as m
        left join lateral (select * from vhs_1708.atlas as s
        where prim=1 and mergedclass=-1
        and jerrorbit=0 and kerrorbit=0
        and q3c_join(m.ra, m.dec,s.ra,s.dec,%0.5f/3600)
        order by q3c_dist(m.ra,m.dec,s.ra,s.dec) asc limit 1)
        as tt on  true  order by xid """ % dist_max,
                               'mytable', (ra, dec, np.arange(
                                   len(dec))), ('ra', 'dec', 'xid'), **kwargs)
    data3 = sqlutil.local_join("""
        select japercormag4,japercormag4_err,
        kapercormag4,kapercormag4_err,ebv from mytable as m
        left join lateral (select * from vhs_1708.gps as s
        where prim=1 and mergedclass=-1
        and jerrorbit=0 and kerrorbit=0
        and q3c_join(m.ra, m.dec,s.ra,s.dec,%0.5f/3600)
        order by q3c_dist(m.ra,m.dec,s.ra,s.dec) asc limit 1)
        as tt on  true  order by xid """ % dist_max,
                               'mytable', (ra, dec, np.arange(
                                   len(dec))), ('ra', 'dec', 'xid'), **kwargs)
    df1 = pd.DataFrame(columns=data1.keys())
    for c in data1.keys():
        df1[c] = data1[c]
    df2 = pd.DataFrame(columns=data2.keys())
    for c in data1.keys():
        df2[c] = data2[c]
    df3 = pd.DataFrame(columns=data3.keys())
    for c in data3.keys():
        df3[c] = data3[c]
    fltr = (df1['japercormag4'] != df1['japercormag4'])
    for c in df2.columns:
        df1.loc[fltr, c] = df2[c][fltr]
    fltr = (df1['japercormag4'] != df1['japercormag4'])
    for c in df3.columns:
        df1.loc[fltr, c] = df3[c][fltr]
    return df1
