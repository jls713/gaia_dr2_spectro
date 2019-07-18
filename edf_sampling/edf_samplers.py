# Script for running MCMC sampling of EDF
# Constructs mock sample, then re-samples the velocities

import sys, os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/py/')
import edf_sampling
import json
import numpy as np
import emcee
from itertools import izip
import pandas as pd
import subprocess
import pickle

# default parameters
names = ['age', 'RcP', 'mass', 'vR', 'vphi', 'vz', 'I', 'l', 'b', 's',
         'Z', 'Teff', 'logg', 'R', 'phi', 'z', 'vlos', 'pm_l', 'pm_b',
         'RA', 'DEC', 'pm_ra', 'pm_dec', 'I', 'JR', 'Lz', 'Jz', 'Rc', 'logl']

def LogL_Zsample(*args):
    return edf_sampling.LogL_sample(*args)

def deccut_fn(l,b,deccut):
    ad = edf_sampling.check_dec(l,b)
    if(deccut<np.pi):
        if(deccut>0.):
            if(ad<deccut-.5*np.pi):
                return 0
        if(deccut<0.):
            if(ad>deccut+.5*np.pi):
                return 0
    return 0

def run_mcmc_sample(
                 Nsamples, mag_band, maglimits,
                 # Isochrone set -- currently probably only Padova will work in all cases
                 which_iso='Padova',
                 # Extra magnitudes to calculate
                 extra_mags=None,
                 # Introduce a colour cut
                 color=None, colorlimits=None, 
                 # Coordinate cuts (in radians) -- remove everything |b|<modbcut
                 # If deccut>0 cut all dec<deccut-pi/2
                 # If deccut<0 cut all dec>deccut+pi/2
                 modbcut=0.,deccut=2.*np.pi,
                 # Spectroscopic parameter cut -- log10Teff
                 logglimits=None,Tefflimits=None,fehlimits=None,
                 # Flags for different options 
                 extinct=True, with_halo=True, interp=False,dered=False,
                 # Debugging and outputting
                 output_file=None,messages=True, extra_magnitudes=None,
                 # Sampler parameters
                 nwalkers=1000,Nburn=None,threads=1,thin=20,pt=False,ntemp=1,asampler=2.,
                 # Adjust magnitude selection with a taper (1-break_grad(mag-break_pos))
                 #break_pos=1000.,break_grad=1000.,
                 # If cutting on error-convolved teff, logg and color -- currently won't do anything
                 tgcolor_errs=np.array([0.,0.,0.])):

    ''' 
	Draws a Monte Carlo sample from an EDF model
	
	Draws Nsamples across the sky between magnitudes maglimits for band mag_band
    '''
    nsamples = int(Nsamples/nwalkers)+1
   
    # Turn on/off halo 
    if(with_halo):
        edf_sampling.turn_on_halo()
    else:
        edf_sampling.turn_off_halo()

    # -- Set arguments
    check = lambda x: np.array([-10000.,10000.]) if x is None else x
    checkcolor = lambda x: np.array(['J','K']) if x is None else color # as a default
    args =[which_iso,mag_band,maglimits,extinct,interp,modbcut,deccut,check(color),
           check(colorlimits),check(Tefflimits),check(logglimits),check(fehlimits)]
    
    # -- Set sampler
    ndim, nwalker = 9, nwalkers
    sampler = emcee.EnsembleSampler(nwalker, ndim,
                                    LogL_Zsample, threads=threads,
                                    args=args)
    if(pt):
        # Parallel tempering
        sampler = emcee.PTSampler(ntemp, nwalker, ndim,los_magbox_LogL_sample_py, 
                                  logp, loglargs=llos_args, threads=threads)
    # Initialize walkers
    # Sample uniform in age, mass
    # Sample gaussian in Z and velocities
    # Sample exp in (unextincted) magnitude
    lomag,himag=maglimits[0],maglimits[1]
    Zlim=[-3.,1.]
    if fehlimits is not None:
        Zlim=fehcut
    if not with_halo:
        Zlim[0]=edf_sampling.minZ()
        Zlim[1]=edf_sampling.maxZ()
    p0 = np.random.uniform(
        low=[0.2, Zlim[0], 0.5, -150., 50., -100., lomag, 0., modbcut],
        high=[12., Zlim[1], 3.5, 150., 350., 100., himag, 2.*np.pi, np.pi /2.],
        size=(nwalker, ndim))
    p0[:,-1]*=(1.-2.*(np.random.uniform(low=0,high=1,size=nwalker)>0.5))
    if fehlimits is not None:
        meanZ,sigZ = 0.,.4
        p0.T[1]=np.random.normal(size=nwalker)*sigZ+meanZ
    meanV,sigV = 0.,50.
    for i in range(3,6):
        p0.T[i]=np.random.normal(size=nwalker)*sigV+meanV
    p0.T[4]+=200.
    p0.T[6]=trunc_exp_rv(lomag,himag,1.,size=len(p0))

    # -- Only choose physically allowed masses for selected age, Z
    for i in np.arange(len(p0)):
        maxmass,minmass=0.,0.
        maxmass = edf_sampling.check_highmass_Z(p0[i][0], p0[i][1],which_iso)
        minmass = edf_sampling.check_lowmass_Z(p0[i][0], p0[i][1],which_iso)
        while(maxmass < 0. or minmass<0. or 
	      edf_sampling.check_radius_positive(p0[i][0], p0[i][1]) == 0 or 
              deccut_fn(p0[i][-2],p0[i][-1],deccut)):
            p0[i][-3] = np.random.uniform(low=lomag, high=himag)
            p0[i][1] = np.random.uniform(low=Zlim[0], high=Zlim[1])
            p0[i][0] = np.random.uniform(low=0., high=12.)
            maxmass = edf_sampling.check_highmass_Z(p0[i][0], p0[i][1],which_iso)
            minmass = edf_sampling.check_lowmass_Z(p0[i][0], p0[i][1],which_iso)
        MAXMASS=3.
        if(maxmass>MAXMASS):
            maxmass=MAXMASS
        p0[i][2] = np.random.uniform(low=minmass, high=maxmass)
        # Check other cuts satisfied
        if(check(colorlimits)[0]>-10. or check(logglimits)[0]>-10. or check(Tefflimits)[0]>-10.):
            while(edf_sampling.check_color_logg_cut(p0[i][0], p0[i][1], p0[i][2], 
                                                 0.,180.,checkcolor(color),
						 check(colorlimits),check(logglimits),check(Tefflimits),
						 which_iso,False)==0):
                    p0[i][2] = np.random.uniform(low=minmass, high=maxmass)
   
    # Now wiggle the walkers a bit 
    err = [0.5,0.1,0.1,10.,10.,10.,0.1,0.01,0.01]
    for i in np.arange(len(p0)):
        extinctMag = edf_sampling.get_extinct(p0[i][-2],p0[i][-1],3.,mag_band)
        if(extinctMag>err[-3]):
            err[-3]=extinctMag
        pp = edf_sampling.LogL_sample(p0[i],*args)
        n=0
        while(np.isinf(pp) and n<10000):
            p0[i]=p0[np.random.randint(len(p0))]+np.random.normal(size=ndim)*err
            pp = edf_sampling.LogL_sample(p0[i],*args)
            n+=1
        if(n==10000):
            print "Can't find start point:",p0[i],deccut_fn(p0[i][-2],p0[i][-1],deccut)

    if(messages): print 'Initial points sampled'
    if(pt):
        p0=np.reshape(p0,(ntemp,nwalker,ndim))

    if(Nburn==None):
        Nburn=2*thin*nsamples
    # -- Run a burn-in
    pos, prob, state = sampler.run_mcmc(p0, Nburn, storechain=False)
    if messages: print 'Number of logl=-inf = '+str(len(prob[np.isinf(prob)]))
    if messages: print("Burnt")
    sampler.reset()

    # -- Sample with thinning and calculating dependent variables
    pos, prob, state = sampler.run_mcmc(pos, nsamples*thin, thin=thin)
    if messages: print ("Sampled")
    flatchain = sampler.flatchain
    lnprob = sampler.lnprobability
    if(pt):
        flatchain=flatchain[0]
        lnprob=lnprob[0]
    extras = np.array(
        map(lambda i: edf_sampling.get_extra_data(i,mag_band,which_iso,False,extinct,interp), sampler.flatchain))
    actions = np.array([edf_sampling.get_actions(np.concatenate((b[4:7], a[3:6])))
                        for a, b in izip(flatchain, extras)])

    nameslist = np.copy(names)
    nameslist[1]="Z"
    nameslist[10]="RcP"
    nameslist[23]=mag_band
    nameslist[6]=mag_band+"0"

    if messages: print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))
    everything = np.vstack((flatchain.T, extras.T,
                     actions.T,lnprob.flatten())).T
    if color is not None:
        if extra_magnitudes is not None:
            extra_magnitudes=np.unique(np.concatenate((color,extra_magnitudes,np.array([mag_band]))))
        else:
            extra_magnitudes=np.unique(np.concatenate((color,np.array([mag_band]))))
        extra_magnitudes = extra_magnitudes[extra_magnitudes!=mag_band]
    
    if extra_magnitudes is not None:
        extra_mags = np.array(
        map(lambda i: edf_sampling.get_extra_magnitudes(i,
                       mag_band, extra_magnitudes,which_iso, False, extinct,interp), flatchain))
        nameslist=np.concatenate((nameslist,extra_magnitudes))
        print extra_magnitudes
        nameslist=np.concatenate((nameslist,np.array([e+'0' for e in list(extra_magnitudes)])))
        everything = np.vstack((everything.T, extra_mags.T)).T
    
    df = pd.DataFrame(everything,columns=nameslist)
    df = df.sample(Nsamples,replace=False).reset_index(drop=True)
    if(output_file):
        df.to_csv(output_file)
    return df

def los_magbox_LogL_sample_py(*args):
    return edf_sampling.los_magbox_LogL_sample(*args)

import scipy.stats

def trunc_exp_rv(lo, hi, scale, size):
    low = 0.
    high=hi-lo
    rnd_cdf = np.random.uniform(scipy.stats.expon.cdf(x=low, scale=scale),
                                scipy.stats.expon.cdf(x=high, scale=scale),
                                size=size)
    return hi-scipy.stats.expon.ppf(q=rnd_cdf, scale=scale)

def logp(p):
    return 1.

def run_mcmc_los(Nsamples, l, b, mag_band, maglimits,
                 # Field radius (in radians)
                 fieldradius=-10.,
                 # Isochrone set -- currently probably only Padova will work in all cases
                 which_iso='Padova',
                 # Extra magnitudes to calculate
                 extra_mags=None,
                 # Introduce a colour cut
                 color=None, colorlimits=None, 
                 # Spectroscopic parameter cut
                 logglimits=None,Tefflimits=None,
                 # Flags for different options 
                 extinct=True, RcPorZ=False, with_halo=True, interp=False,dered=False,
                 # Debugging and outputting
                 output_file=None,messages=True, extra_magnitudes=None,
                 # Sampler parameters
                 nwalkers=1000,Nburn=None,threads=1,thin=20,pt=False,ntemp=1,asampler=2.,
                 # Adjust magnitude selection with a taper (1-break_grad(mag-break_pos))
                 break_pos=1000.,break_grad=1000.,
                 # If cutting on error-convolved teff, logg and color
                 tgcolor_errs=np.array([0.,0.,0.])):

    ''' 
	Draws a Monte Carlo sample from an EDF model
	
	Draws Nsamples in a field of radius 
        fieldradius (in radians) at (l,b) (in radians) between magnitudes maglimits for
        band mag_band
    '''
    nsamples = int(Nsamples/nwalkers)+1

    if(with_halo):
        edf_sampling.turn_on_halo()
    else:
        edf_sampling.turn_off_halo()

    if(np.fabs(b)<fieldradius/4.):
        b+=fieldradius/2.

    # -- Initialise walkers
    ndim, nwalker = 7 + 2*(fieldradius>0), nwalkers
    # Parameters for log-likelihood
    check = lambda x: np.array([-10000.,10000.]) if x is None else x
    checkcolor = lambda x: np.array(['J','K']) if x is None else color
    llos_args =[np.array([l,b]), mag_band, maglimits, which_iso, RcPorZ, 
                checkcolor(color), check(colorlimits), check(logglimits), check(Tefflimits), 
                np.array([1*extinct,1*dered]), fieldradius,interp,
                np.array([break_pos,break_grad]),tgcolor_errs]
    ## Samplers 
    sampler = emcee.EnsembleSampler(nwalker, ndim,los_magbox_LogL_sample_py, 
                                    args=llos_args, threads=threads,a=asampler)
    if(pt):
        # Parallel tempering
        sampler = emcee.PTSampler(ntemp, nwalker, ndim,los_magbox_LogL_sample_py, 
                                  logp, loglargs=llos_args, threads=threads)
    # Either sample in birth radius or metallicity -- birth radius sampling obsolete but good for checks
    RorZlim=[0.,18.]
    if(RcPorZ==False):
        RorZlim=[-1.,0.6]
        if not with_halo:
            RorZlim[0]=edf_sampling.minZ()
            RorZlim[1]=edf_sampling.maxZ()
    
    # Sample uniform in age, mass
    # Sample gaussian in Z and velocities
    # Sample exp in (unextincted) magnitude
    meanZ,sigZ = 0.,.4
    meanV,sigV = 0.,50.
    #     age  Z           mass  vr     vphi    vz   mag
    lo = [0.5, RorZlim[0], 0.5, -150., -200., -150., maglimits[0]]
    hi = [12., RorZlim[1], 2.,   150.,  200.,  150., maglimits[1]]
    if(fieldradius>0.):
        # also sample in l,b 
        hi=np.concatenate((hi,[l,b]))
        lo=np.concatenate((lo,[l,b]))
    # Initial sample
    p0 = np.random.uniform(low=lo,high=hi,size=(nwalker*ntemp, ndim))
    p0.T[1]=np.random.normal(size=(nwalker*ntemp))*sigZ+meanZ # metallicity
    for i in range(3,6):
        # velocities
        p0.T[i]=np.random.normal(size=(nwalker*ntemp))*sigV+meanV
    p0.T[4]+=200. # shift vphi
    # magnitudes
    p0.T[6]=trunc_exp_rv(maglimits[0],maglimits[1],1.,size=len(p0))

    # -- Only choose physically allowed masses for selected age, Z
    for i in np.arange(len(p0)):
        maxmass,minmass=0.,0.
        check_highmass_fn = edf_sampling.check_highmass_Z
        check_lowmass_fn = edf_sampling.check_lowmass_Z
        if(RcPorZ==True):
            check_highmass_fn = edf_sampling.check_highmass
            check_lowmass_fn = edf_sampling.check_lowmass
        maxmass = check_highmass_fn(p0[i][0], p0[i][1],which_iso)
        minmass = check_lowmass_fn(p0[i][0], p0[i][1],which_iso)

        while(maxmass < 0. or minmass<0. or 
              edf_sampling.check_radius_positive(p0[i][0], p0[i][1]) == 0):
            p0[i][1] = np.random.uniform(low=RorZlim[0], high=RorZlim[1])
            maxmass = check_highmass_fn(p0[i][0], p0[i][1],which_iso)
            minmass = check_lowmass_fn(p0[i][0], p0[i][1],which_iso)
        MAXMASS=3.
        if(maxmass>MAXMASS):
            maxmass=MAXMASS
        p0[i][2] = np.random.uniform(low=minmass, high=maxmass)
        # Check other cuts satisfied
        if(check(colorlimits)[0]>-10. or check(logglimits)[0]>-10. or check(Tefflimits)[0]>-10.):
            while(edf_sampling.check_color_logg_cut(p0[i][0], p0[i][1], p0[i][2], 
                                                 0.,180.,checkcolor(color),
						 check(colorlimits),check(logglimits),check(Tefflimits),
						 which_iso,False)==0):
                p0[i][2] = np.random.uniform(low=minmass, high=maxmass)

    # Got reasonable starting points but now want to shuffle slightly to ensure
    err = [0.5,0.1,0.1,10.,10.,10.,0.1]
    # Shuffle by extinction if large
    max_distance=3.
    extinctBand = edf_sampling.get_extinct(l,b,max_distance,mag_band)
    if(extinct>err[-1]):
        err[-1]=extinctBand
    if(fieldradius>0.):
        err=np.concatenate((err,[0.,0.]))
    # Shuffle
    for i in np.arange(len(p0)):
        pp = los_magbox_LogL_sample_py(p0[i],*llos_args)
        n,maxn=0,10000
        while(np.isinf(pp) and n<maxn):
            p0[i]=p0[np.random.randint(len(p0))]+np.random.normal(size=ndim)*err
            pp = los_magbox_LogL_sample_py(p0[i],*llos_args)
            n+=1
        if(n==maxn):
            print 'No finite LogL found:', p0[i]
    ## Now sample l and b if needed
    if(fieldradius>0.):
        rad = np.random.uniform(size=len(p0))*fieldradius
        thet = np.random.uniform(size=len(p0))*2.*np.pi
        p0.T[7]=l+np.cos(thet)*rad/np.cos(b)
        p0.T[8]=b+np.sin(thet)*rad

    if(messages): print 'Initial points sampled'
    if(pt):
        p0=np.reshape(p0,(ntemp,nwalker,ndim))
    # -- Run a burn-in
    Nb=16000
    if(Nburn):
        Nb=Nburn
    else:
        Nb=nsamples*thin*2
    pos, prob, state = sampler.run_mcmc(p0, Nb, storechain=False)
    if(messages):
        print 'Number of logl=-inf = '+str(len(prob[np.isinf(prob)]))
        print("Burnt")
    sampler.reset()

    # -- Sample with thinning and calculating dependent variables
    pos, prob, state = sampler.run_mcmc(pos, nsamples*thin, thin=thin)
    flatchain = sampler.flatchain
    lnprob = sampler.lnprobability
    if(pt):
        flatchain=flatchain[0]
        lnprob=lnprob[0]

    extras = np.array(
        map(lambda i: edf_sampling.get_extra_data(i if fieldradius>0. else np.append(i,[l,b]), 
                       mag_band, which_iso, RcPorZ, extinct,interp), flatchain))
    actions = np.array([edf_sampling.get_actions(np.concatenate((j[4:7], i[3:6])))
                        for i,j in izip(flatchain, extras)])

    nameslist = np.copy(names)
    if(RcPorZ==False):
        nameslist[1]="Z"
        nameslist[10]="RcP"
    nameslist[23]=mag_band
    nameslist[6]=mag_band+"0"
    
    LB = np.reshape(np.tile([l,b],len(flatchain)),(len(flatchain),2))
    everything = np.vstack((flatchain.T, LB.T, extras.T,
                 actions.T,lnprob.flatten())).T
    if fieldradius>0.:
        everything = np.vstack((flatchain.T, extras.T,
                 actions.T,lnprob.flatten())).T

    if color is not None:
        if extra_magnitudes is not None:
            extra_magnitudes=np.unique(np.concatenate((color,extra_magnitudes,np.array([mag_band]))))
        else:
            extra_magnitudes=np.unique(np.concatenate((color,np.array([mag_band]))))
        extra_magnitudes = extra_magnitudes[extra_magnitudes!=mag_band]
    
    if extra_magnitudes is not None:
        extra_mags = np.array(
        map(lambda i: edf_sampling.get_extra_magnitudes(i if fieldradius>0. else np.append(i,[l,b]),
                       mag_band, extra_magnitudes,which_iso, RcPorZ, extinct,interp), flatchain))
        nameslist=np.concatenate((nameslist,extra_magnitudes))
        print extra_magnitudes
        nameslist=np.concatenate((nameslist,np.array([e+'0' for e in list(extra_magnitudes)])))
        everything = np.vstack((everything.T, extra_mags.T)).T

    if(messages): print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))
    df=pd.DataFrame(everything,columns=nameslist)
    df = df.sample(Nsamples,replace=False).reset_index(drop=True)
    if(output_file):
        df.to_csv(output_file)
    return df

from sklearn.neighbors import KDTree

def add_errors_to_sample(sample,data,
                         data_tree_fields = ['Teff_K','logg_K','Met_N_K_DR5'],
                         sample_tree_fields = ['Teff','logg','Z'],
                         err_fields = {'Teff':'eTeff_K',
                                       'logg':'elogg_K',
                                       'Z':'eMet_N_K',
                                       'vlos':'eHRV'},
                         leaf_size=10):
    '''
        Add errors to a mock sample
        sample -- sample catalogue
        data -- data
        data_tree_fields -- the data entries used to construct the tree
        sample_tree_fields -- the sample entries corresponding to the data used
                              to generate the tree
        err_fields -- a dictionary of the fields to add errors to in the keys
                      and the corresponding error fields in the data
    '''
    # Build tree
    data_tr = data[data_tree_fields]
    m,s = data_tr.mean(), data_tr.std()
    data_tr = (data_tr-m)/s
    kdt = KDTree(data_tr, leaf_size=leaf_size, metric='euclidean')
    RAVE_dd = (sample[sample_tree_fields].values-m.values)/s.values

    index = kdt.query(RAVE_dd, k=1, return_distance=False)
    for k,e in err_fields.items():
        sample['e'+k]=data[e].take(index.flatten()).values
        sample[k+'_e']=np.random.normal(size=len(sample))*sample['e'+k]+sample[k]

    return sample

def add_errors_to_sample_astrometry_covariance(sample,data,
                         data_tree_fields = ['BVmag','Vmag','l','b'],
                         sample_tree_fields = ['BV','V','l','b'],
                         err_fields = {'B':'Bmag','V':'Vmag'},
                         leaf_size=10):
    '''
        Add errors to a mock sample
        sample -- sample catalogue
        data -- data
        data_tree_fields -- the data entries used to construct the tree
        sample_tree_fields -- the sample entries corresponding to the data used
                              to generate the tree
        err_fields -- a dictionary of the fields to add errors to in the keys
                      and the corresponding error fields in the data
        Will also use tree fields to assign covariant astrometry errors
    '''
    # Build tree
    data_tr = data[data_tree_fields]
    m,s = data_tr.mean(), data_tr.std()
    data_tr = (data_tr-m)/s
    kdt = KDTree(data_tr, leaf_size=leaf_size, metric='euclidean')
    sample_dd = (sample[sample_tree_fields].values-m.values)/s.values

    sample['parallax'] = 1./sample['s']

    index = kdt.query(sample_dd, k=1, return_distance=False)
    ## First photometry
    for k,e in err_fields.items():
        sample['e'+k]=data[e].take(index.flatten()).values
        sample[k+'_e']=np.random.normal(size=len(sample))*sample['e'+k]+sample[k]
    ## Now astrometry
    cov = np.array([[data.parallax_error.values**2,
        data.parallax_pmra_corr.values*data.parallax_error.values*data.pmra_error.values,
        data.parallax_pmdec_corr.values*data.parallax_error.values*data.pmdec_error.values],
       [data.parallax_pmra_corr.values*data.parallax_error.values*data.pmra_error.values,
        data.pmra_error.values**2,
        data.pmra_pmdec_corr.values*data.pmra_error.values*data.pmdec_error.values],
       [data.parallax_pmdec_corr.values*data.parallax_error.values*data.pmdec_error.values,
        data.pmra_pmdec_corr.values*data.pmra_error.values*data.pmdec_error.values,
        data.pmdec_error.values**2]]).T
    #DD = np.array(map(lambda x: np.random.multivariate_normal(np.zeros(3),x),
    #        cov[index.flatten()]))
    DD = np.einsum('jik,jk->ji',np.linalg.cholesky(cov[index.flatten()]),
                                np.random.normal(size=(len(sample),3)))
    for i,p in enumerate(['parallax','pm_ra','pm_dec']):
        sample[p+'_e']=DD[:,i]+sample[p]
        sample['e'+p]=np.sqrt(cov[index.flatten()][:,i,i])
    sample['parallax_pmra_corr']=cov[index.flatten()][:,1,0]/np.sqrt(cov[index.flatten()][:,0,0]*cov[index.flatten()][:,1,1])
    sample['parallax_pmdec_corr']=cov[index.flatten()][:,2,0]/np.sqrt(cov[index.flatten()][:,0,0]*cov[index.flatten()][:,2,2])
    sample['pmra_pmdec_corr']=cov[index.flatten()][:,1,2]/np.sqrt(cov[index.flatten()][:,2,2]*cov[index.flatten()][:,1,1])
    return sample

def compute_velocities_actions_with_errors(sample):
    '''
        compute the velocities and actions of a sample with erroneous
        quantities.
        sample must have fields RA, DEC, s_e, vlos_e, pm_ra_e, pm_dec_e

        NEED TO RUN edf_sampling.setup first

    '''

    fields =['PMl_e','PMb_e','R_e','phi_e','z_e','vR_e','vphi_e','vz_e','JR_e','Lz_e','Jz_e','Rc_e']
    for f in fields:
        sample[f]=np.ones(len(sample))*-9999.

    for i in range(len(sample)):
        if(sample['s_e'][i]>0. and sample['s_e'][i]==sample['s_e'][i]):
            X = edf_sampling.process_data(
                np.array([sample['RA'][i],     sample['DEC'][i],
                          sample['s_e'][i],    sample['vlos_e'][i],
                          sample['pm_ra_e'][i],sample['pm_dec_e'][i]]))
            sample['PMl_e'][i],sample['PMb_e'][i]=X[4],X[5]
            sample['R_e'][i],sample['phi_e'][i],sample['z_e'][i]=X[6],X[7],X[8]
            sample['vR_e'][i],sample['vphi_e'][i],sample['vz_e'][i]=X[9],X[10],X[11]
            sample['JR_e'][i],sample['Lz_e'][i],sample['Jz_e'][i],sample['Rc'][i]=X[12],X[13],X[14],X[15]

    return sample


def compute_velocities_with_errors(sample):
    '''
        compute the velocities and actions of a sample with erroneous
        quantities.
        sample must have fields RA, DEC, s_e, vlos_e, pm_ra_e, pm_dec_e
    '''

    fields =['PMl_e','PMb_e','R_e','phi_e','z_e','vR_e','vphi_e','vz_e','JR_e','Lz_e','Jz_e','Rc_e']
    for f in fields:
        sample[f]=np.ones(len(sample))*-9999.

    for i in range(len(sample)):
        if(sample['s_e'][i]>0. and sample['s_e'][i]==sample['s_e'][i]):
            X = edf_sampling.process_data(
                np.array([sample['RA'][i],     sample['DEC'][i],
                          sample['s_e'][i],    sample['vlos_e'][i],
                          sample['pm_ra_e'][i],sample['pm_dec_e'][i]]))
            sample['PMl_e'][i],sample['PMb_e'][i]=X[4],X[5]
            sample['R_e'][i],sample['phi_e'][i],sample['z_e'][i]=X[6],X[7],X[8]
            sample['vR_e'][i],sample['vphi_e'][i],sample['vz_e'][i]=X[9],X[10],X[11]
            sample['JR_e'][i],sample['Lz_e'][i],sample['Jz_e'][i],sample['Rc'][i]=X[12],X[13],X[14],X[15]

    return sample
