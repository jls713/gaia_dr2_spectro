from subprocess import call
import pandas as pd
from astroquery.vizier import Vizier
import numpy as np


def download_passbands():
    '''
        Downloads the quantum efficiency for SDSS, PanSTARRS,
        Landolt B&V, WISE, 2MASS
        and theoretical Gaia G, G_RP, G_BP, G_RVS
    '''
    #filters = ['SDSSu', 'SDSSg', 'SDSSr', 'SDSSi', 'SDSSz',
               # 'LandoltB', 'LandoltV',
    #           'WISEW1', 'WISEW2',
    #           '2MASSJ', '2MASSH', '2MASSKs']
    #for f in filters:
    #    call(['rm', f])
    #    call(['wget', '-nH', '-nd', '-np',
    #          'http://www.mso.anu.edu.au/~brad/Filters/' + f])

    for f in ['u','g','r','i','z']:
        call(['rm', 'SDSS'+f])
        call(['wget', '-nH', '-nd', '-np',
	      'http://classic.sdss.org/dr3/instruments/imager/filters/%s.dat'%f,
	       '-O', 'SDSS%s'%f])
        sdss_bands = pd.read_csv('SDSS%s'%f, sep=r'\s+',skiprows=6,
				names=['Angstrom','pass','extended','exxtended0','airmass1'])
        sdss_bands[['Angstrom', 'pass']].to_csv('SDSS%s' % f, sep=' ',
                                           header=['0.', ''], index=False)
    for f in ['W1', 'W2']:
        call(['rm', 'WISE'+f])
        call(['wget', '-nH', '-nd', '-np',
	      'http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/figures/RSR-%s.txt'%f,
	       '-O', 'WISE%s'%f])
        wise_bands = pd.read_csv('WISE%s'%f, sep=r'\s+',skiprows=2,
				names=['Angstrom','pass','uncertainty'])
        wise_bands['Angstrom']*=10000.
        wise_bands[['Angstrom', 'pass']].to_csv('WISE%s' % f, sep=' ',
                                           header=['0.', ''], index=False)

    for i,f in enumerate(['J', 'H', 'Ks']):
        call(['rm', '2MASS'+f])
        call(['wget', '-nH', '-nd', '-np',
              'http://www.ipac.caltech.edu/2mass/releases/second/doc/sec3_1b1.tbl%i.html'%(i+12),
	       '-O', '2MASS%s'%f])
        tmass_bands = pd.read_csv('2MASS%s'%f, sep=r'\s+',skiprows=11, skipfooter=7,
				names=['count','Angstrom','pass'])
	## Need to transform to photon counting
        tmass_bands['pass']/=tmass_bands['Angstrom']
        tmass_bands['pass']/=np.max(tmass_bands['pass'])
        tmass_bands['Angstrom']*=10000.
        tmass_bands[['Angstrom', 'pass']].to_csv('2MASS%s' % f, sep=' ',
                                           header=['0.', ''], index=False)
    return

    ubv = Vizier.get_catalogs('J/AJ/131/1184/table3')[0]
    for fltr in ['U', 'B', 'V']:
        tbl = ubv[['%slambda' % fltr, '%sph' % fltr]].to_pandas()
        tbl.to_csv('Landolt%s' % fltr, sep=' ',
                   header=['0.', ''], index=False)

    call(['rm', 'PSg', 'PSr', 'PSi', 'PSz', 'PSy'])
    call(['wget', '-nH', '-nd', '-np',
          'http://ipp.ifa.hawaii.edu/ps1.filters/apj425122t3_mrt.txt',
          '-O', 'PanSTARRS.txt'])
    ps_bands = pd.read_csv('PanSTARRS.txt', sep=r'\s+', skiprows=26,
                           names=['nm', 'whole', 'g', 'r', 'i', 'z',
                                  'y', 'w', 'aero', 'ray', 'mol'])
    ps_bands['Angstrom'] = ps_bands['nm'] * 10.
    for fld in ['g', 'r', 'i', 'z', 'y']:
        ps_bands[['Angstrom', fld]].to_csv('PS%s' % fld, sep=' ',
                                           header=['0.', ''], index=False)

    # Gaia pre-launch
    call(['rm', 'GaiaG_pre', 'GaiaBP_pre', 'GaiaRP_pre', 'GaiaRVS_pre'])
    call(['wget', '-nH', '-nd', '-np',
          'https://www.cosmos.esa.int/documents/29201/302420/'
          'normalisedPassbands.txt/a65b04bd-4060-44fa-be36-91975f2bd58a',
          '-O', 'GaiaBands.txt'])

    gaia_bands = pd.read_csv('GaiaBands.txt', sep=r'\s+', skiprows=[1])
    gaia_bands['Angstrom'] = gaia_bands['LAMBDA'] * 10.
    for fld in ['G', 'BP', 'RP', 'RVS']:
        gaia_bands[['Angstrom', fld]].to_csv('Gaia%s_pre' % fld, sep=' ',
                                             header=['0.', ''], index=False)

    # Gaia DR2 processing
    call(['rm', 'GaiaG_dr2', 'GaiaBP_dr2', 'GaiaRP_dr2'])
    call(['wget', '-nH', '-nd', '-np',
          'https://www.cosmos.esa.int/documents/29201/1645651/'
          'GaiaDR2_Passbands_ZeroPoints.zip/'
          '49cdce41-8eee-655d-7ed2-4e7a83598c1d',
          '-O', 'GaiaDR2_Passbands_ZeroPoints.zip'])
    call(['unzip','GaiaDR2_Passbands_ZeroPoints.zip'])
    call(['rm','GaiaDR2_Passbands_ZeroPoints.zip'])

    gaia_bands = pd.read_csv('GaiaDR2_Passbands.dat', sep=r'\s+',
                             names=['LAMBDA','G','Gerr',
                             'BP','BPerr','RP','RPerr'])
    gaia_bands['Angstrom'] = gaia_bands['LAMBDA'] * 10.
    for fld in ['G', 'BP', 'RP']:
    	gaia_bands.loc[gaia_bands[fld]>1.,fld]=0.
        gaia_bands[['Angstrom', fld]].to_csv('Gaia%s_dr2' % fld, sep=' ',
                                             header=['0.', ''], index=False)

    # Gaia revised
    call(['rm', 'GaiaG_revised', 'GaiaBP_revised', 'GaiaRP_revised'])
    call(['wget', '-nH', '-nd', '-np',
          'https://www.cosmos.esa.int/documents/29201/1645651/'
          'GaiaDR2_Revised_Passbands_ZeroPoints.zip/'
          '54db454f-69cb-ea0c-15be-f1b1f597f191',
          '-O', 'GaiaDR2_Revised_Passbands_ZeroPoints.zip'])
    call(['unzip','GaiaDR2_Revised_Passbands_ZeroPoints.zip'])
    call(['rm','GaiaDR2_Revised_Passbands_ZeroPoints.zip'])

    gaia_bands = pd.read_csv('GaiaDR2_RevisedPassbands.dat', sep=r'\s+',
                             names=['LAMBDA','G','Gerr',
                             'BP','BPerr','RP','RPerr'])
    gaia_bands['Angstrom'] = gaia_bands['LAMBDA'] * 10.
    for fld in ['G', 'BP', 'RP']:
    	gaia_bands.loc[gaia_bands[fld]>1.,fld]=0.
        gaia_bands[['Angstrom', fld]].to_csv('Gaia%s_revised' % fld, sep=' ',
                                             header=['0.', ''], index=False)


if __name__ == '__main__':
    download_passbands()
