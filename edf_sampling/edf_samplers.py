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

names = ['age', 'RcP', 'mass', 'vR', 'vphi', 'vz', 'I', 'l', 'b', 's',
         'Z', 'Teff', 'logg', 'R', 'phi', 'z', 'vlos', 'pm_l', 'pm_b',
         'RA', 'DEC', 'pm_ra', 'pm_dec', 'J-K', 'J', 'K', 'I', 'JR', 'Lz', 'Jz', 'Rc', 'logl']

names_Z = ['age', 'Z', 'mass', 'vR', 'vphi', 'vz', 'I', 'l', 'b', 's',
           'RcP', 'Teff', 'logg', 'R', 'phi', 'z', 'vlos', 'pm_l', 'pm_b',
           'RA', 'DEC', 'pm_ra', 'pm_dec', 'JminusK', 'J','K','I0','JR', 'Lz', 'Jz', 'Rc',
           'logl']

names_resample = ['age', 'vR', 'vphi', 'vz', 'RcP']

names_V = ['age', 'Z', 'mass', 'vR', 'vphi', 'vz', 'V']
names_VRCP = ['age', 'RcP', 'mass', 'vR', 'vphi', 'vz', 'V0','l','b', 's',
              'Z', 'Teff', 'logg', 'R', 'phi', 'z', 'vlos', 'pm_l', 'pm_b',
              'RA', 'DEC', 'pm_ra', 'pm_dec','JminusK','V','JR', 'Lz', 'Jz', 'Rc',
           'logl']

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

def run_mcmc_sample(Nsamples,Nburn=None,which='BaSTI',thin=20,output_file=None,band="I",maglimits=[8.,13.],extinct=True,with_halo=False,add_APASS=False,interp=False,modbcut=0.,deccut=2.*np.pi,JKcut=np.array([0.,0.]),Teffcut=np.array([-1e6,1e6]),loggcut=np.array([-1e6,1e6]),fehcut=np.array([-1e6,1e6])):
    '''
        maglimits -- mag limits for initial sample
        Cut in modb in radians
        deccut>0 we cut everything with dec<deccut-PI/2
        deccut<0 we cut everything with dec>deccut+PI/2
        Teffcut we keep everything with Teffcut[0]<log10(Teff)<Teffcut[1]
        loggcut we keep everything with loggcut[0]<logg<loggcut[1]
    '''
    # with open('config.json') as data_file:
    #     data = json.load(data_file)
    #     which = str(data['isochrones'])

    Zlim=[-3.,1.]
    if not with_halo:
        Zlim[0]=edf_sampling.minZ()
        Zlim[1]=edf_sampling.maxZ()

    # -- Initialise walkers
    ndim, nwalker = 9, 5000
    args =[which,band,extinct,interp,modbcut,deccut,JKcut,Teffcut,loggcut,fehcut]
    sampler = emcee.EnsembleSampler(nwalker, ndim,
                                    LogL_Zsample, threads=8,
                                    args=args)
    lomag,himag=maglimits[0],maglimits[1]
    p0 = np.random.uniform(
        low=[0.2, Zlim[0], 0.5, -150., 50., -100., lomag, 0., -np.pi / 2.],
        high=[12., Zlim[1], 3.5, 150., 350., 100., himag, 2.*np.pi, np.pi /2.],
        size=(nwalker, ndim))
    meanZ,sigZ = 0.,.4
    meanV,sigV = 0.,50.
    p0.T[1]=np.random.normal(size=nwalker)*sigZ+meanZ
    for i in range(3,6):
        p0.T[i]=np.random.normal(size=nwalker)*sigV+meanV
    p0.T[4]+=200.
    p0.T[6]=trunc_exp_rv(lomag,himag,1.,size=len(p0))

    MAXMASS=3.

    # -- Only choose physically allowed masses for selected age, Z
    for i in np.arange(len(p0)):
        maxmass,minmass=0.,0.
        maxmass = edf_sampling.check_highmass_Z(p0[i][0], p0[i][1],which)
        minmass = edf_sampling.check_lowmass_Z(p0[i][0], p0[i][1],which)

        while(maxmass < 0. or minmass<0. or edf_sampling.check_radius_positive(p0[i][0], p0[i][1]) == 0 or deccut_fn(p0[i][-2],p0[i][-1],deccut)):
            p0[i][2] = np.random.uniform(low=0., high=MAXMASS)
            p0[i][-3] = np.random.uniform(low=lomag, high=himag)
            p0[i][1] = np.random.uniform(low=Zlim[0], high=Zlim[1])
            p0[i][0] = np.random.uniform(low=0., high=12.)
            maxmass = edf_sampling.check_highmass_Z(p0[i][0], p0[i][1],which)
            minmass = edf_sampling.check_lowmass_Z(p0[i][0], p0[i][1],which)

    for i in np.arange(len(p0)):
        err = [0.5,0.1,0.1,10.,10.,10.,0.1,0.01,0.01]
        extinctMag = edf_sampling.get_extinct(p0[i][-2],p0[i][-1],3.,band)
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

    print("Found suitable initial conditions")
    for i in p0:
        if(np.isinf(edf_sampling.LogL_sample(i,*args))):
            print "-inf log-likelihood for params: ", i

    if(Nburn==None):
        Nburn=10*Nsamples
    # -- Run a burn-in
    pos, prob, state = sampler.run_mcmc(p0, Nburn, storechain=False)
    print 'Number of logl=-inf = '+str(len(prob[np.isinf(prob)]))
    print("Burnt")
    sampler.reset()

    # -- Sample with thinning and calculating dependent variables
    pos, prob, state = sampler.run_mcmc(pos, Nsamples, thin=thin)
    print ("Sampled")
    flatchain = sampler.flatchain
    lnprob = sampler.lnprobability
    extras = np.array(
        map(lambda i: edf_sampling.get_extra_data(i,band,which,False,extinct,interp), sampler.flatchain))
    if(add_APASS):
        apass = np.array(
            map(lambda i: edf_sampling.get_APASS(i,band,which,False,extinct,interp), sampler.flatchain))
    extras=np.vstack((extras.T,apass.T)).T
    actions = np.array([edf_sampling.get_actions(np.concatenate((b[4:7], a[3:6])))
                        for a, b in izip(flatchain, extras)])

    nameslist = names[:]
    nameslist[1]="Z"
    nameslist[10]="RcP"
    mag_band=band
    nameslist[26]=mag_band
    nameslist[6]=mag_band+"0"
    if(mag_band=='J'):
	nameslist[24]='J_extra'
    if(mag_band=='K'):
	nameslist[25]='K_extra'
    if(add_APASS):
        nameslist.insert(27,'B')
        nameslist.insert(28,'V')
    print nameslist
    if(output_file):
        np.savetxt(output_file, np.vstack((flatchain.T, extras.T, actions.T,
                               lnprob)).T,
               header="\t".join([p for p in nameslist]),
               delimiter='\t', comments='')

    print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))
    everything = np.vstack((flatchain.T, extras.T,
                     actions.T,lnprob.flatten())).T
    return pd.DataFrame(everything,columns=nameslist)

def los_magbox_LogL_sample_py(*args):
    return los_magbox_LogL_sample(*args)

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

def run_mcmc_los(Nsamples, nwalkers, l, b, which_iso, lomag=-10., himag=15., mag_band="I",JminusKcut=np.array([-100.,100.]), loggcut=np.array([-2000.,2000.]),Teffcut=np.array([.01,20000.]),extinct=True, RcPorZ=False, Nburn=None, with_halo=True, output_file=None,messages=True,threads=1,thin=20,pt=False,ntemp=1,fieldradius=-10.,interp=False,asampler=2.,break_pos=1000.,break_grad=1000.,dered=False,tgjk_errs=np.array([0.,0.,0.])):

    # with open('config.json') as data_file:
    #     data = json.load(data_file)
    #     which_iso = str(data['isochrones'])

    if(with_halo):
        edf_sampling.turn_on_halo()
    else:
        edf_sampling.turn_off_halo()

    if(np.fabs(b)<fieldradius/4.):
        b+=fieldradius/2.

    # -- Initialise walkers
    ndim, nwalker = 7, nwalkers
    if(fieldradius>0.):
        ndim+=2
    llos_args =[l, b, mag_band, np.array([lomag, himag]), which_iso, RcPorZ, JminusKcut, loggcut, Teffcut, np.array([1*extinct,1*dered]), fieldradius,interp,np.array([break_pos,break_grad]),tgjk_errs]
    sampler = emcee.EnsembleSampler(nwalker, ndim,los_magbox_LogL_sample_py, args=llos_args, threads=threads,a=asampler)
    if(pt):
        sampler = emcee.PTSampler(ntemp, nwalker, ndim,los_magbox_LogL_sample_py, logp, loglargs=llos_args, threads=threads)
    RorZlim=[0.,18.]
    if(RcPorZ==False):
        RorZlim=[-1.,0.6]
        if not with_halo:
            RorZlim[0]=edf_sampling.minZ()
            RorZlim[1]=edf_sampling.maxZ()
    MAXMASS=3.
    # Sample uniform in age, mass
    # Sample gaussian in Z and velocities
    # Sample exp in magnitude
    lo = [0.5, RorZlim[0], 0.5, -150., -200., -150., lomag]
    hi  = [12., RorZlim[1], 2., 150., 200., 150., himag]
    if(fieldradius>0.):
        hi=np.concatenate((hi,[l,b]))
        lo=np.concatenate((lo,[l,b]))
    p0 = np.random.uniform(low=lo,high=hi,size=(nwalker*ntemp, ndim))
    meanZ,sigZ = 0.,.4
    meanV,sigV = 0.,50.
    p0.T[1]=np.random.normal(size=(nwalker*ntemp))*sigZ+meanZ
    for i in range(3,6):
        p0.T[i]=np.random.normal(size=(nwalker*ntemp))*sigV+meanV
    p0.T[4]+=200.
    p0.T[6]=trunc_exp_rv(lomag,himag,1.,size=len(p0))

    # -- Only choose physically allowed masses for selected age, Z
    for i in np.arange(len(p0)):
        maxmass,minmass=0.,0.
        if(RcPorZ==False):
            maxmass = edf_sampling.check_highmass_Z(p0[i][0], p0[i][1],which_iso)
            minmass = edf_sampling.check_lowmass_Z(p0[i][0], p0[i][1],which_iso)
        else:
            maxmass = edf_sampling.check_highmass(p0[i][0], p0[i][1],which_iso)
            minmass = edf_sampling.check_lowmass(p0[i][0], p0[i][1],which_iso)
        while(maxmass < 0. or minmass<0. or edf_sampling.check_radius_positive(p0[i][0], p0[i][1]) == 0):
            p0[i][1] = np.random.uniform(low=RorZlim[0], high=RorZlim[1])
            if(RcPorZ==False):
                maxmass = edf_sampling.check_highmass_Z(p0[i][0], p0[i][1],which_iso)
                minmass = edf_sampling.check_lowmass_Z(p0[i][0], p0[i][1],which_iso)
            else:
                maxmass = edf_sampling.check_highmass(p0[i][0], p0[i][1],which_iso)
                minmass = edf_sampling.check_lowmass(p0[i][0], p0[i][1],which_iso)
        if(maxmass>MAXMASS):
            maxmass=MAXMASS
        p0[i][2] = np.random.uniform(low=minmass, high=maxmass)
        if(JminusKcut[0]>-10. or loggcut[0]>-10. or Teffcut[0]>-10.):
            while(edf_sampling.check_JK_logg_cut(p0[i][0], p0[i][1], p0[i][2], 0.,180.,JminusKcut,loggcut,Teffcut,which_iso,False)==0):
                p0[i][2] = np.random.uniform(low=minmass, high=maxmass)
    llarray=np.zeros(len(p0))
    for i in np.arange(len(p0)):
        llarray[i] = los_magbox_LogL_sample_py(p0[i],*llos_args)
    maxll = np.max(llarray)-60.
    err = [0.5,0.1,0.1,10.,10.,10.,0.1]
    extinctH = edf_sampling.get_extinct_H(l,b,3.)
    if(extinctH>err[-1]):
        err[-1]=extinctH
    if(fieldradius>0.):
        err=np.concatenate((err,[0.,0.]))
    for i in np.arange(len(p0)):
        pp = los_magbox_LogL_sample_py(p0[i],*llos_args)
        n=0
        while(np.isinf(pp) and n<10000):
            p0[i]=p0[np.random.randint(len(p0))]+np.random.normal(size=ndim)*err
            pp = los_magbox_LogL_sample_py(p0[i],*llos_args)
            n+=1
        if(n==10000):
            print p0[i]
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
        Nb=Nsamples*thin*2
    pos, prob, state = sampler.run_mcmc(p0, Nb, storechain=False)
    if(messages):
        print 'Number of logl=-inf = '+str(len(prob[np.isinf(prob)]))
        print("Burnt")
    sampler.reset()

    # -- Sample with thinning and calculating dependent variables
    pos, prob, state = sampler.run_mcmc(pos, Nsamples*thin, thin=thin)
    flatchain = sampler.flatchain
    lnprob = sampler.lnprobability
    if(pt):
        flatchain=flatchain[0]
        lnprob=lnprob[0]

    extras = np.array(
        map(lambda i: edf_sampling.get_extra_data(i if fieldradius>0. else np.append(i,[l,b]), mag_band, which_iso, RcPorZ, extinct,interp), flatchain))
    actions = np.array([edf_sampling.get_actions(np.concatenate((j[4:7], i[3:6])))
                        for i,j in izip(flatchain, extras)])

    nameslist = names
    if(RcPorZ==False):
        nameslist[1]="Z"
        nameslist[10]="RcP"
    nameslist[26]=mag_band
    nameslist[6]=mag_band+"0"

    LB = np.reshape(np.tile([l,b],len(flatchain)),(len(flatchain),2))
    everything = np.vstack((flatchain.T, LB.T, extras.T,
                 actions.T,lnprob.flatten())).T
    if fieldradius>0.:
        everything = np.vstack((flatchain.T, extras.T,
                 actions.T,lnprob.flatten())).T

    if(output_file):
        np.savetxt(output_file, everything, header="\t".join([p for p in nameslist]),
               delimiter='\t', comments='')

    if(messages): print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))
    return pd.DataFrame(everything,columns=nameslist)


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

    sample.parallax = 1./sample.s

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
    mean = np.array([data.parallax.values,data.pmra.values,data.pmdec.values]).T
    DD = np.array(map(lambda x: np.random.multivariate_normal(x[0],x[1]),
            izip(mean[index.flatten()],cov[index.flatten()]))).T
    for i,p in enumerate(['parallax','pm_ra','pm_dec']):
        sample[p+'_e']=DD[i]
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
