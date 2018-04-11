import numpy as np
import isodist_js as isodist
import emcee
from numpy.linalg import inv
#==============================================================================

def init(types="All",wemap=True):
    isodist.init_isochrone(types,1,0.5)
    if(wemap):
        isodist.load_emap()
    isodist.load_prior('2018')

class distance_data:
    def __init__(self,mag,mag_err,magstr,inputs,cov,prior,which,
                 parallax=[0.,-1.],ak=0.):
        if(len(cov)==6):
            s=3
            alpha=0
        else:
            s=4
            alpha=1
        covm = np.zeros((s,s))
        covm[0][0]=cov[0]
        covm[0][1]=cov[1]
        covm[1][0]=cov[1]
        covm[0][2]=cov[2]
        covm[2][0]=cov[2]
        covm[1][1]=cov[3+alpha]
        covm[1][2]=cov[4+alpha]
        covm[2][1]=cov[4+alpha]
        covm[2][2]=cov[5+alpha*2]
        if(alpha):
            covm[0][3]=cov[3]
            covm[3][0]=cov[3]
            covm[1][3]=cov[6]
            covm[3][1]=cov[6]
            covm[2][3]=cov[8]
            covm[3][2]=cov[8]
            covm[3][3]=cov[9]
        icovm = inv(covm)
        self.icov = np.array([icovm[0][0],icovm[0][1],icovm[0][2],icovm[1][1],icovm[1][2],icovm[2][2]])
        if(alpha):
            self.icov=np.insert(self.icov,3,icovm[0][3])
            self.icov=np.insert(self.icov,-1,icovm[1][3])
            self.icov=np.append(self.icov,icovm[2][3])
            self.icov=np.append(self.icov,icovm[3][3])
        self.cov_alphafree = np.array([covm[0][0],covm[0][1],covm[0][2],covm[1][1],covm[1][2],covm[2][2]])
        self.mag=mag
        self.mag_err=mag_err
        self.magstr=magstr
        self.inputs=inputs
        self.prior=prior
        self.which=which
        self.cov=cov
        self.ak=ak
        self.parallax = parallax[0]
        self.parallax_error = parallax[1]

import time

def distance_emcee(data,extinct=False,alpha=False,nwalkers=30,
                   use_extinct=False,nsteps=1000,asampler=1.3):
    '''
        If extinct, marginalize over extinction.
        If alpha, marginalize over alpha
    '''
    def logl(x,data):
        if(x[2]<0. or x[1]<0. or x[1]>13.5 or x[0]<-3.5 or x[0]>1.):
            return -np.inf
        ddd = isodist.distance_pdf(x[-1],x[:-1],data.mag,data.inputs,data.mag_err,data.icov,data.prior,data.which,data.magstr)
        if(np.isnan(ddd)):
            return -np.inf
        else:
            return ddd-(extinct)*0.5*(data.ak-np.exp(x[3])/8.93)**2/0.05**2
    start = time.time()
    ddata = isodist.prob_distance(data.mag,data.inputs,data.mag_err,
                                  data.cov_alphafree,data.prior,
                                  data.which,data.magstr,use_extinct,
                                  0.,False,data.parallax,data.parallax_error)
    print ddata
    print 'Time=',time.time()-start
    nnum=4+extinct+alpha
    sampler = emcee.EnsembleSampler(nwalkers,nnum,logl,args=[data],a=asampler)
    feh_err = np.sqrt(data.cov[0])
    if(feh_err!=feh_err or np.isnan(feh_err)):
        feh_err=0.05
    p = [ddata[11],ddata[7],ddata[9],ddata[1]]
    sp= [feh_err,ddata[8],ddata[10],ddata[2]]

    if(extinct):
        p.insert(3,np.log(data.ak*8.93))
        sp.insert(3,0.05)
    if(alpha):
        p.insert(-1,data.inputs[-1])
        sp.insert(-1,1./np.sqrt(data.icov[-1]))
    print p,sp


    start = time.time()

    pos, prob, state = sampler.run_mcmc([np.random.normal(loc=p,scale=sp) for i in range(nwalkers)],nsteps)
    pos00 = pos[(prob!=-np.inf)&(prob==prob)]

    pos=pos00[np.random.randint(0,pos00.shape[0],len(pos))]
    sampler.reset()
    # pos, prob, state = sampler.run_mcmc(pos,nsteps)
    # pos00 = pos[prob>np.max(prob)-100.]
    # pos=pos00[np.random.randint(0,pos00.shape[0],len(pos))]
    # sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos,nsteps)
    pos00 = pos[prob>np.max(prob)-20.]
    pos=pos00[np.random.randint(0,pos00.shape[0],len(pos))]
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos,nsteps)
    sampler.reset()
    sampler.run_mcmc(pos,nsteps)
    print 'Time=',time.time()-start
    print np.mean(sampler.acceptance_fraction)
    tefflogg = [isodist.get_teff_logg(i,j,k,data.which,True,l) for i,j,k,l in zip(sampler.flatchain.T[0],sampler.flatchain.T[1],sampler.flatchain.T[2],sampler.flatchain.T[4])]
    return sampler.flatchain,tefflogg,sampler.flatlnprobability

def distance_emcee_parallel_tempering(data,extinct=False,alpha=False,nwalkers=30,ntemps=20,use_extinct=False,nsteps=1000,asampler=2.):
    '''
        If extinct, marginalize over extinction.
        If alpha, marginalize over alpha
    '''
    def logp(x):
        if(x[2]<0. or x[1]<0. or x[1]>13.5 or x[0]<-3.5 or x[0]>1.):
            return -np.inf
        else:
            return 0.
    def logl(x,data):
        ddd = isodist.distance_pdf(x[-1],x[:-1],data.mag,data.inputs,data.mag_err,data.icov,data.prior,data.which,data.magstr)
        if(np.isnan(ddd)):
            return -np.inf
        else:
            return ddd-(extinct)*0.5*(data.ak-np.exp(x[3])/8.93)**2/0.05**2

    ddata = isodist.prob_distance(data.mag,data.inputs,data.mag_err,data.cov_alphafree,data.prior,data.which,data.magstr,use_extinct,0.,False,data.parallax,data.parallax_error)
    print ddata
    nnum=4+extinct+alpha
    sampler = emcee.PTSampler(ntemps,nwalkers,nnum,logl,logp,loglargs=[data],a=asampler)

    feh_err = np.sqrt(data.cov[0])
    if(feh_err!=feh_err or np.isnan(feh_err)):
        feh_err=0.05
    p = [ddata[11],ddata[7],ddata[9],ddata[1]]
    sp= [feh_err,ddata[8],ddata[10],ddata[2]]
    if(extinct):
        p.insert(3,np.log(data.ak*8.93))
        sp.insert(3,0.05)
    if(alpha):
        p.insert(-1,data.inputs[-1])
        sp.insert(-1,1./np.sqrt(data.icov[-1]))

    print p,sp
    p0=[[np.random.normal(loc=p,scale=sp) for i in range(nwalkers)] for j in range(ntemps)]
    pos, prob, like = sampler.run_mcmc(p0,nsteps)
    pos00 = pos[(prob!=-np.inf)&(prob==prob)]
    pos=np.array([pos00[np.random.randint(0,pos00.shape[0],len(pos[0]))] for i in range(len(pos))])
    sampler.reset()
    pos, prob, like = sampler.run_mcmc(pos,nsteps)
    sampler.reset()
    pos, prob, like = sampler.run_mcmc(pos,nsteps*10,thin=10)
    print np.mean(sampler.acceptance_fraction)
    # np.max(sampler.acor)
    print np.shape(sampler.chain)
    print np.shape(sampler.flatchain[0,...])
    tefflogg = [isodist.get_teff_logg(i,j,k,data.which,True,l) for i,j,k,l in zip(sampler.flatchain[0,...].T[0],sampler.flatchain[0,...].T[1],sampler.flatchain[0,...].T[2],sampler.flatchain[0,...].T[4])]
    return sampler.flatchain[0,...],tefflogg#,sampler.lnprobability.reshape(ntemps,nwalkers*nsteps)[0]

def distance_emcee_es(data,extinct=False,alpha=False,nwalkers=30,use_extinct=False,nsteps=1000,asampler=1.3):
    '''
        If extinct, marginalize over extinction.
        If alpha, marginalize over alpha
    '''

    def logl(x, data):
        if(x[2] < 0. or x[2]>1. or x[1]<0. or x[1]>13.5 or x[0]<-3.5 or x[0]>1.):
            return -np.inf
        ddd = isodist.distance_pdf_es(x[-1],x[:-1],data.mag,data.inputs,data.mag_err,data.icov,data.prior,data.which,data.magstr)
        if(np.isnan(ddd)):
            return -np.inf
        else:
            return ddd-(extinct)*0.5*(data.ak-np.exp(x[3])/8.93)**2/0.05**2
    ddata = isodist.prob_distance(data.mag,data.inputs,data.mag_err,data.cov_alphafree,data.prior,data.which,data.magstr,use_extinct,0.,False,data.parallax,data.parallax_error)

    print ddata
    nnum=4+extinct+alpha
    sampler = emcee.EnsembleSampler(nwalkers,nnum,logl,args=[data],a=asampler)
    feh_err = np.sqrt(data.cov[0])
    if(feh_err!=feh_err or np.isnan(feh_err)):
        feh_err=0.05
    es = isodist.get_es(ddata[11],ddata[7],ddata[9],data.which)
    p = [ddata[11],ddata[7],es,ddata[1]]
    sp= [feh_err,2.,0.02,ddata[2]]
    if(extinct):
        p.insert(3,np.log(data.ak*8.93))
        sp.insert(3,0.05)
    if(alpha):
        p.insert(-1,data.inputs[-1])
        sp.insert(-1,1./np.sqrt(data.icov[-1]))
    print p,sp



    pos, prob, state = sampler.run_mcmc([np.random.normal(loc=p,scale=sp) for i in range(nwalkers)],nsteps)
    pos00 = pos[prob!=-np.inf]
    pos=pos00[np.random.randint(0,pos00.shape[0],len(pos))]
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos,nsteps)
    pos00 = pos[prob>np.max(prob)-100.]
    pos=pos00[np.random.randint(0,pos00.shape[0],len(pos))]
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos,nsteps)
    pos00 = pos[prob>np.max(prob)-20.]
    pos=pos00[np.random.randint(0,pos00.shape[0],len(pos))]
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos,nsteps)
    sampler.reset()
    sampler.run_mcmc(pos,2000)
    print np.mean(sampler.acceptance_fraction)
    tefflogg = [isodist.get_teff_logg(i,j,k,data.which,True) for i,j,k in zip(sampler.flatchain.T[0],sampler.flatchain.T[1],sampler.flatchain.T[2])]
    return sampler.flatchain,tefflogg,sampler.flatlnprobability

def test():
    init(types="Padova", wemap=True)

    mags = np.array([12.1239996, 13.4119997, 11.59599972])
    errmags = np.array([0.029, 0.026, 0.021])
    mag_str = np.array(['H', 'J', 'K'])

    data=np.array([-0.133,3.65738258,2.14995,1.57964417,0.01856532])#,0.])
    cov = np.array([0.1*0.1,0.,0.,0.1*0.1,0.,0.1*0.1])
    prior=True
    w="Padova"

    ddata = distance_data(mags,errmags,mag_str,data,cov,prior,w)
    # ddd = distance_emcee_parallel_tempering(ddata,extinct=True,use_extinct=True,alpha=True,nsteps=500,asampler=2.,nwalkers=100,ntemps=10)
    # ddd = distance_emcee(ddata,extinct=True,use_extinct=True,alpha=True,nsteps=500,asampler=2.,nwalkers=500)
    ddd = distance_emcee_es(ddata,extinct=False,use_extinct=False,
                            alpha=False,nsteps=500,asampler=2.,nwalkers=200)
    return ddd
    print np.mean(ddd,axis=0),np.std(ddd,axis=0)
    ddd = distance_emcee(ddata,extinct=True)
    print np.mean(ddd,axis=0),np.std(ddd,axis=0)
    cov = np.array([ERRZ*ERRZ,0.,0.,0.,ERRTEFF*ERRTEFF,0.,0.,ERRLOGG*ERRLOGG,0.,0.05])
    data = np.array([np.float64(Z),np.float64(TEFF),np.float64(LOGG),np.float64(l),np.float64(b),0.])
    ddata = distance_data(mags,errmags,mag_str,data,cov,prior,w)
    ddd = distance_emcee(ddata,extinct=True,alpha=True)
    print np.mean(ddd,axis=0),np.std(ddd,axis=0)

    # mags = np.array([14.2257,14.7692,14.0935])
    # Z = -0.04
    # TEFF = 3.65741
    # LOGG = 2.67457
    # # ERRZ=0.02
    # # ERRTEFF=0.005
    # # ERRLOGG=0.06
    # l=0.43570463
    # b=0.84198128
    # errmags=np.array([0.022,0.027,0.021])
    # prior=True
    # w="Padova"
    # mag_str=np.array(["H","J","K"])
    # data = np.array([np.float64(Z),np.float64(TEFF),np.float64(LOGG),np.float64(l),np.float64(b)])
    # data_errs = np.array([np.float64(ERRZ),np.float64(ERRTEFF),np.float64(ERRLOGG)])
    # print isodist.prob_distance(mags,data,errmags,data_errs,prior,w,mag_str,False,0.,False)
    # cov = np.array([0.01,7.96797003e-06,3.93751893e-05,3.51248734e-05,5.34178871e-05,3.62309873e-03])


    return ddd

def test2():
    init(types="Dartmouth",wemap=False)
    mags = np.array([10.83918594 , 11.34971984  ,10.75608826])
    mags = np.array([11.10099983  ,11.77200031  ,10.92500019])
    Z = -0.041
    TEFF = 3.67416432
    LOGG = 2.61385
    l=3.33877442
    b=0.06705567
    errmags=np.array([0.03,0.027,0.019])
    prior=True
    w="Dartmouth"
    mag_str=np.array(["H","J","K"])
    data = np.array([np.float64(Z),np.float64(TEFF),np.float64(LOGG),np.float64(l),np.float64(b)])
    # data_errs = np.array([np.float64(ERRZ),np.float64(ERRTEFF),np.float64(ERRLOGG)])
    # print isodist.prob_distance(mags,data,errmags,data_errs,prior,w,mag_str,False,0.,False)
    cov = np.array([1.43135510e-05 ,  2.18613244e-05 ,  9.62819429e-05 ,  8.42088351e-05,
   3.19634694e-04 ,  5.64131993e-03])
    ddata = distance_data(mags,errmags,mag_str,data,cov,prior,w,0.169)
    ddd = distance_emcee(ddata,extinct=True)
    return ddd

if __name__ == '__main__':
    test()
