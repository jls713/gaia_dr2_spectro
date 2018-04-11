##=============================================================================
import numpy as np
import fitsio
import pandas as pd
import matplotlib.pyplot as plt
import h5py
import seaborn as sns
from scipy.interpolate import RectBivariateSpline
import sys
from numpy.lib.recfunctions import append_fields
##=============================================================================

def andy_OK(data):
    non_fucked = (data.TEFF_ASPCAP>3500.)&(data.TEFF_ASPCAP<5500.)&(data.LOGG_ASPCAP>0.)&(data.LOGG_ASPCAP<3.9)
    data = data[non_fucked].reset_index(drop=True)
    return data

def read_fits_to_pandas(data_file):
    output_data=pd.DataFrame()
    input_data = fitsio.FITS(data_file)[1].read().byteswap().newbyteorder()
    for i in input_data.dtype.descr:
        if(isinstance(input_data[i[0]][0],np.ndarray)):
            for c,j in enumerate(input_data[i[0]].T):
                if(isinstance(j,np.ndarray)):
                    continue
                else:
                    output_data[i[0]+str(c)]=j
        else:
            output_data[i[0]] = input_data[i[0]]
    return output_data

data_file = '/data/arc/research/apogee/regularized-results/tc-cse-regularized-apogee-catalog.fits.gz'
apogee_data = read_fits_to_pandas(data_file)
apogee_data = andy_OK(apogee_data)

##=============================================================================

Alpha_Elements = ['O','Ne','Mg','Si','S','Ar','Ca','Ti']
SolarRadius=8.3
SunAge=4.
GalaxyAge = 13.7

class chem_evo_data:
    ## Class to handle the data and plotting of the results of the chemical
    ## evolution modelling
    def __init__(self,filename):
        ''' Load in the results of the chemical evolution model from hdf5 '''
        data = h5py.File(filename,'r')
        self.R = data['R'][()][0].T
        self.t = data['t'][()][0].T
        self.Mgas = np.dstack(data['Mgas'][()])[0]
        self.Mstar = np.dstack(data['Mstar'][()])[0]
        self.Z = np.dstack(data['Z'][()])[0]
        self.SFR = data['SFR'][()][0].T
        self.Inflow = data['Inflow'][()][0].T
        self.SNIa = data['SNIa'][()][0].T
        self.SNII = data['SNII'][()][0].T
        self.elements = [i for i in data.keys() if i not in ['R','t','Mgas','Z','parameters','SFR','Inflow','SNIa','SNII','Mstar']]
        for e in self.elements:
            if e in Alpha_Elements:
                self.elements += [u'alpha']
                print self.elements
                break
        self.abund = np.recarray((len(self.R),len(self.t)),
                                  dtype=[(i.encode('ascii','ignore'),np.float)
                                    for i in self.elements])
        Nalpha=0
        for e in self.elements[:-1]:
            self.abund[e]=np.dstack(data[e][()])[0].T
            if e in Alpha_Elements:
                if Nalpha==0:
                    self.abund['alpha']=self.abund[e]
                else:
                    self.abund['alpha']+=self.abund[e]
                Nalpha+=1
        if(Nalpha>0):
            self.abund['alpha']/=float(Nalpha)
        data.close()

    def plot_radial(self,el,el2=None,time=GalaxyAge,color='k',show_gradient=False):
        ''' Plot the radial distribution of the element el at time time. if
            el2 != None el-el2 is plotted '''
        dat = self.abund[el]
        if(el2):
            dat=dat-self.abund[el2]
        rbs = RectBivariateSpline(self.R,self.t,dat)
        R = np.linspace(self.R[0],self.R[-1],100)
        a = rbs(R,time)
        plt.plot(R,a,color=color)
        plt.xlabel(r'$R/\mathrm{kpc}$')
        if(el2):
            plt.ylabel(r'$[\mathrm{%s}/\mathrm{%s}]$'%(el,el2))
        else:
            plt.ylabel(r'$[\mathrm{%s}]$'%(el))
        if(show_gradient):
            grad=.5*(rbs(SolarRadius+1.,GalaxyAge)-rbs(SolarRadius-1.,GalaxyAge))
            plt.annotate(r'$%0.3f/\mathrm{dex\,kpc}$'%(grad),xy=(0.95,0.95),xycoords='axes fraction',ha='right',va='top')
        plt.tight_layout()

    def plot_radius_range(self,el,timerange=np.arange(0.1,12.1,1.),el2=None):
        ''' Plots the radial distribution of element el for a range of times
            timerange. if el2 != None el-el2 is plotted '''
        for t in timerange:
            self.plot_radial(el,el2,time=t)
        self.plot_radial(el,el2,time=GalaxyAge-SunAge,color=sns.color_palette()[2],show_gradient=True)
        plt.plot(SolarRadius,0.,'.',ms=10)

    def plot_time(self,el,el2=None,radius=SolarRadius,color='k'):
        ''' Plot the distribution in time of the element el at radius radius. If
            el2 != None el-el2 is plotted '''
        dat = self.abund[el]
        if(el2):
            dat=dat-self.abund[el2]
        rbs = RectBivariateSpline(self.R,self.t,dat)
        t = np.linspace(self.t[0],self.t[-1],100)
        a = rbs(radius,t)
        plt.plot(t,a.T,color=color)
        plt.xlabel(r'$t/\mathrm{Gyr}$')
        if(el2):
            plt.ylabel(r'$[\mathrm{%s}/\mathrm{%s}]$'%(el,el2))
        else:
            plt.ylabel(r'$[\mathrm{%s}]$'%(el))
        plt.tight_layout()

    def plot_time_range(self,el,radiusrange=np.arange(1.,15.,1.),el2=None):
        ''' Plots the distribution in time of element el for a range of radii
            radiusrange. if el2 != None el-el2 is plotted '''
        for r in radiusrange:
            self.plot_time(el,el2,radius=r)
        self.plot_time(el,el2,radius=SolarRadius,color=sns.color_palette()[2])
        plt.plot(self.t[-1]-SunAge,0.,'.',ms=10)

    def plot_elements_against_time(self,el,el2,el_u=None,el_u2=None,radius=SolarRadius,color='k'):
        ''' Plot element el against el2 for all times at fixed radius. If el_u
            is not None plot el-el_u and similarly if el_u2 is not None plot
            el2-el_u2 '''
        dat = self.abund[el]
        if(el_u):
            dat=dat-self.abund[el_u]
        rbs = RectBivariateSpline(self.R,self.t,dat)
        t = np.linspace(self.t[0],self.t[-1],100)
        a = rbs(radius,t)
        dat = self.abund[el2]
        if(el_u2):
            dat=dat-self.abund[el_u2]
        rbs = RectBivariateSpline(self.R,self.t,dat)
        b = rbs(radius,t)
        plt.plot(a[0],b[0],color=color)
        if(el_u2):
            plt.ylabel(r'$[\mathrm{%s}/\mathrm{%s}]$'%(el2,el_u2))
        else:
            plt.ylabel(r'$[\mathrm{%s}]$'%(el2))
        if(el_u):
            plt.xlabel(r'$[\mathrm{%s}/\mathrm{%s}]$'%(el,el_u))
        else:
            plt.xlabel(r'$[\mathrm{%s}]$'%(el))
        plt.tight_layout()

    def plot_element_range(self,el,el2,el_u=None,el_u2=None,radiusrange=np.arange(1.,15.,1.)):
        ''' Plot element el against element el2 for all times at a range of
            radii radiusrange. If el_u is not None plot el-el_u and similarly
            if el_u2 is not None plot el2-el_u2 '''
        for r in radiusrange:
            self.plot_elements_against_time(el,el2,el_u,el_u2,radius=r)
        self.plot_elements_against_time(el,el2,el_u,el_u2,radius=SolarRadius,color=sns.color_palette()[2])
        plt.plot(0.,0.,'.',ms=10)

    def plot_abundance_matrix(self):
        ''' Plot the APOGEE data for each element with the chemical evolution
            tracks overplotted '''
        elts = [i for i in self.elements if i not in ['H','He']]
        f=plt.figure(figsize=[10.,10.])
        n_plts=len(elts)
        kk=1
        for e in elts:
            plt.subplot(6,n_plts,kk)
            self.plot_time_range(e,el2='H')
            kk+=1
        for e in elts:
            plt.subplot(6,n_plts,kk)
            self.plot_radius_range(e,el2='H')
            kk+=1
        for e in elts:
            plt.subplot(6,n_plts,kk)
            self.plot_time_range(e,el2='Fe')
            kk+=1
        for e in elts:
            plt.subplot(6,n_plts,kk)
            self.plot_radius_range(e,el2='Fe')
            kk+=1
        for e in elts:
            plt.subplot(6,n_plts,kk)
            plt.plot(apogee_data.FE_H,apogee_data[e.upper()+'_H']-apogee_data.FE_H,'.',alpha=0.1,color=sns.color_palette()[0])
            self.plot_element_range('Fe',e,'H','Fe')
            kk+=1
        for e in elts:
            plt.subplot(6,n_plts,kk)
            plt.plot(apogee_data.O_H-apogee_data.FE_H,apogee_data[e.upper()+'_H']-apogee_data.FE_H,'.',alpha=0.1,color=sns.color_palette()[0])
            self.plot_element_range('O',e,'Fe','Fe')
            kk+=1
        plt.tight_layout()

    def summary_plot(self):
        '''
            Plot a summary plot for the chemical evolution model.
            The star formation rate, inflow rate, SNIa and SNII rates at the
            solar radius are plotted along with observational constraints from
            Guesten & Mezger et al. (1982) (for the SFR)
            Portinari et al. (1998) (for the inflow)
            and Li et al. (2011) (for the SN rates)
        '''
        plt.figure(figsize=[3.,4.])
        plt.subplot(2,1,1)
        plt.plot(self.t,self.SFR,color='k',label='SFR')
        l,=plt.plot(self.t,self.Inflow,color='k',ls='dashed',label='Inflow')
        l.set_dashes((2,1))
        plt.legend(loc='upper left', bbox_to_anchor=(0.1, 1.0))
        plt.xlabel(r'$t/\,\mathrm{Gyr}$')
        plt.ylabel(r'Rate /$\,\mathrm{M}_\odot\,\mathrm{pc}^{-2}\,\mathrm{Gyr}^{-1}$')
        plt.ylim(0.,np.max(self.SFR)*1.4)

        sfr = [6.,4.] # Guesten & Mezger et al. (1982)
        inf = [0.9,0.6] # Portinari et al. (1998)
        plt.errorbar([11.8],[sfr[0]],yerr=[sfr[1]],color=sns.color_palette()[2],fmt='*',markersize=4)
        plt.errorbar([11.8],[inf[0]],yerr=[inf[1]],color=sns.color_palette()[2],fmt='*',markersize=4)

        ax=plt.twinx()
        plt.plot(self.t,self.SNIa,color='k',ls='dotted',label='SNIa')
        l,=plt.plot(self.t,self.SNII,color='k',ls='dashed',label='SNII')
        l.set_dashes((4,1))
        plt.ylabel(r'Rate /$\,\mathrm{pc}^{-2}\,\mathrm{Gyr}^{-1}$')
        plt.legend(loc='upper right', bbox_to_anchor=(0.9, 1.0))
        plt.ylim(0.,plt.ylim()[1]*1.4)

        typeII = [1.555,0.285] # Li et al. (2011)
        typeIa = [0.54,0.11] # Li et al. (2011)
        ## Assuming 15kpc disc
        conv = 10./np.pi/15.**2
        plt.errorbar([12.2],[typeIa[0]*conv],yerr=[typeIa[1]*conv],color=sns.color_palette()[1],fmt='o',markersize=2)
        plt.errorbar([12.2],[typeII[0]*conv],yerr=[typeII[1]*conv],color=sns.color_palette()[1],fmt='o',markersize=2)
        plt.xlim(0.,self.t[-1]*1.05)
        plt.subplot(2,1,2)
        plt.plot(self.R,self.Mstar[-1],color='k',label='Stars')
        l,=plt.plot(self.R,self.Mgas[-1],color='k',ls='',label='Gas')
        l.set_dashes((2,1))
        plt.legend(loc='upper right', bbox_to_anchor=(0.9, 1.0))
        plt.errorbar([SolarRadius,SolarRadius],[43.,13.],yerr=[5.,3.],fmt='o',markersize=2)
        plt.xlim(2.,16.)
        plt.semilogy()
        plt.xlabel(r'$R/\,\mathrm{kpc}$')
        plt.ylabel(r'Surface Density$/\,\mathrm{M}_\odot\,\mathrm{pc}^{-2}$')
        plt.tight_layout()

    def paint(self,R,t,el):
        ''' Returns the abundances in list el at a range of R and t '''
        results = []
        rbs_h = RectBivariateSpline(self.R,self.t,self.abund['H'])
        for e in el:
            rbs = RectBivariateSpline(self.R,self.t,self.abund[e])
            results+=[rbs(R,t,grid=False)-rbs_h(R,t,grid=False)]
        return results


if __name__=="__main__":
    filename='tmp.h5'
    if(sys.argv[1]):
        filename=sys.argv[1]
    Data = chem_evo_data(filename)
    f=plt.figure()
    Data.summary_plot()
    plt.savefig('summary.pdf')
    f=plt.figure(figsize=[4.,7.])
    plt.subplot('311')
    Data.plot_time_range('Fe',el2='H')
    plt.subplot('312')
    Data.plot_radius_range('Fe',el2='H')
    plt.subplot('313')
    Data.plot_element_range('Fe','alpha','H','Fe')
    plt.savefig('tmp.pdf')
    plt.clf()
    # Data.plot_abundance_matrix()
    # plt.savefig('tmp2.png',dpi=400)
