#include "boost_python_interface.h"
#include "tables_aa.h"
#include "potential.h"
#include "iso_grid.h"
#include "distances.h"
#include "edf.h"
#include "sf.h"
#include "sf_impl.h"
#include "mass_function.h"
INITIALIZE_EASYLOGGINGPP
#include "chemevo/inc/hdf5_reader.h"
using namespace H5;

double dfsf(const VecDoub& v,double *sigma,double Vc){
    return exp(-.5*(pow(v[0]/sigma[0],2)+pow((v[1]-Vc)/sigma[1],2)+pow(v[2]/sigma[2],2)));
}

void printZeros(int n){
    for (int i = 0; i < n; ++i) std::cout<<"0. ";
    std::cout<<"-1e10\n";
}

namespace edf_sampling{
    sb15_edf *EDF;

    std::unique_ptr<isochrone_grid<isochrone_johnson>> iso;
    std::unique_ptr<isochrone_grid<isochrone_padova>> iso_Padova;
    std::unique_ptr<isochrone_grid<isochrone_dartmouth>> iso_Dartmouth;
    std::unique_ptr<DistanceCalculator<isochrone_johnson>> isoD;
    std::unique_ptr<DistanceCalculator<isochrone_padova>> isoD_Padova;
    std::unique_ptr<DistanceCalculator<isochrone_dartmouth>> isoD_Dartmouth;

    std::shared_ptr<extinction_law> EL;
    std::shared_ptr<sfd_extinction_map> EM;

    std::shared_ptr<MassFunction> Mass_Function;

    Selection_function *RSF;
    std::vector<Selection_function*> sf_vec;

    Potential_JS *Pot;

    Actions_AxisymmetricFudge_InterpTables *Tab;
    std::unique_ptr<galaxy_prior> prior;
    int NCA=0;

    rand_uniform *rn;
    rand_gaussian *rnGauss;
    rand_exponential *rnExp;

    VecDoub mc_samples, data_samples;

    void preamble(){
        rn = new rand_uniform(time(0));
        rnGauss = new rand_gaussian(1.,time(0)+1);
        rnExp = new rand_exponential(1.,time(0)+2);
    }

    void load_selection_function(std::string sfn){
        RSF = sf_types[sfn]();
    }
    void load_selection_function_stack(std::vector<std::string> sfn){
        for(auto sf: sfn)
            sf_vec.push_back(sf_types[sf]());
        RSF = new Stacked_Binned_selection_function(sf_vec);
    }
    void init_isochrone(const std::string & which,int thin=1,double feh_err=1e-5,double thin_mag=-1.){
        if(which=="BaSTI"){
            iso = make_unique<isochrone_grid<isochrone_johnson>>(which,thin,feh_err,thin_mag);
            isoD = make_unique<DistanceCalculator<isochrone_johnson>>(&*iso);
        }
        else if(which=="Padova"){
            iso_Padova = make_unique<isochrone_grid<isochrone_padova>>(which,thin,feh_err,thin_mag);
            isoD_Padova = make_unique<DistanceCalculator<isochrone_padova>>(&*iso_Padova);
        }
        else if(which=="Dartmouth"){
            iso_Dartmouth = make_unique<isochrone_grid<isochrone_dartmouth>>(which,thin,feh_err,thin_mag);
            isoD_Dartmouth = make_unique<DistanceCalculator<isochrone_dartmouth>>(&*iso_Dartmouth);
        }
        else if(which=="All"){
            iso = make_unique<isochrone_grid<isochrone_johnson>>("BaSTI",thin,feh_err,thin_mag);
            iso_Padova = make_unique<isochrone_grid<isochrone_padova>>("Padova",thin,feh_err,thin_mag);
            iso_Dartmouth = make_unique<isochrone_grid<isochrone_dartmouth>>("Dartmouth",thin,feh_err,thin_mag);
            isoD = make_unique<DistanceCalculator<isochrone_johnson>>(&*iso);
            isoD_Padova = make_unique<DistanceCalculator<isochrone_padova>>(&*iso_Padova);
            isoD_Dartmouth = make_unique<DistanceCalculator<isochrone_dartmouth>>(&*iso_Dartmouth);
        }
        EL = std::make_shared<schlafly2017_extinction_law>();
    }

    void setup(bool ext_load=true, bool iso_load=true){
        auto params = parameters();
    	preamble();

        std::string acts_file = params["actions"]["dir"];
        std::string potential = params["potential"];
        VecDoub solar_motion = params["solar_motion"];
        bool with_halo = params["with_halo"];
        std::string params_file = params["edf_params"];
        std::string isochrones = params["isochrones"]["type"];

        std::ifstream infile(acts_file);
        bool acts_calculated=false;
        if(infile.good()==true)
            acts_calculated=true;
        infile.close();
        Pot = new GalPot(potential);

        Tab = new Actions_AxisymmetricFudge_InterpTables(Pot,acts_file,!acts_calculated,
            params["actions"]["rmin"],
            params["actions"]["rmax"],
            params["actions"]["nrgrid"],
            params["actions"]["ngrid"]);
    	std::cout<<"InterpTables loaded"<<std::endl;

        EDF = new sb15_edf(Pot,Tab,{0.,0.},solar_motion); // pass nothing params as we read from file below
    	std::cout<<"EDF loaded"<<std::endl;
        EDF->readParams(params_file);
        if(!with_halo)EDF->TurnOffHalo();
        else EDF->TurnOnHalo();

        EL = std::make_shared<schlafly2017_extinction_law>();

        if(ext_load)
            EM = std::make_shared<sfd_extinction_map>(&*EL,solar_motion);
        if(iso_load)
            init_isochrone(isochrones,1,params["isochrones"]["feh_spacing"]);

    	// Load prior
    	auto prior_type = params["prior"];
        if(prior_type=="Binney"){
            prior = make_unique<binney_prior>(solar_motion);
            std::cout<<"Binney 2014 prior loaded for EDF fitting\n";
        }
        else if(prior_type=="2018"){
            prior = make_unique<new_prior_2018>(solar_motion);
            std::cout<<"New 2018 prior loaded for EDF fitting\n";
        }
        else if(prior_type=="2018_broad"){
            prior = make_unique<new_prior_2018_broad_age>(solar_motion);
            std::cout<<"New 2018 broad age prior loaded for EDF fitting\n";
        }
        else{
            prior = make_unique<galaxy_prior>(solar_motion);
            std::cout<<"Flat prior loaded for EDF fitting\n";
        }

        Mass_Function = std::make_shared<MassFunction>(params);
    }
    VecDoub get_params(void){return EDF->returnParams();}
    void    set_params(VecDoub A){return EDF->setParams(A);}
    void    set_params_json(std::string A){return EDF->setParams(A);}
    double get_Vc(double R){return EDF->pot->Vc(R);}

    void load_emap(void){
	EM = std::make_shared<sfd_extinction_map>(&*EL,conv::StandardSolarPAUL);
    }

    VecDoub LogL_sample(double *f,std::string which="BaSTI", std::string band="I", bool extinct=true,bool interp=false, double modbcut=0.,double deccut=PI,VecDoub JKcut={0.,0.}, VecDoub Teffcut = {-1e6,1e6}, VecDoub loggcut = {-1e6,1e6},VecDoub FeHcut={-1e6,1e6}){
        //
        // Calculates log-likelihood given f
        // f = tau, Z, M, vR, vp,  vz, Mag, l, b
        //
        // Cut in modb in radians
        // deccut>0 we cut everything with dec<deccut-PI/2
        // deccut<0 we cut everything with dec>deccut+PI/2
        // JKcut = {Cut value (i.e. keep everything greater than cut),
        //          only apply cut for stars with fabs(b)<JKcut[1]}
        // Teffcut = Keep all Teffcut[0]<log10(Teff)<Teffcut[1]
        // loggcut = Keep all loggcut[0]<logg<loggcut[1]
        // -- Extract quantities from f and check if allowed

        double tau = f[0];
        if(tau>=EDF->get_taum() or tau<0.){
            return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
        }

        double Z = f[1];
        // double RcP = f[1];
        if(Z<-3. or Z>1. or Z<FeHcut[0] or Z>FeHcut[1]) {
            return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
        }

        double M = f[2];
        if(M<0.) {
            return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
        }
        double Mag = f[6];
        double l = f[7];
        double b = f[8];

        if(l<0. or l>2.*PI or b<-PI/2. or b>PI/2. or fabs(b)<modbcut) {
            return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
        }
        if(deccut<PI){
            VecDoub Eq = conv::GalacticToEquatorial({l,b,1.});
            if(deccut>0.){
                if(Eq[1]<deccut-.5*PI) return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
            }
            else if(deccut<0.){
                if(Eq[1]>deccut+.5*PI) return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
            }
        }
        // -- Find nearest isochrone point or interpolate
        double s,JK;
        std::vector<std::string> mag_sf = RSF->mag_band();
        VecDoub sf_mags;
        if(interp){
            VecDoub ite;
            std::vector<std::string> req_bands = {band,"J","K"};
            req_bands.insert(req_bands.end(),mag_sf.begin(),mag_sf.end());
            if(which=="BaSTI")
                ite = iso->interp_es(Z,tau,M,req_bands);
            else if(which=="Padova")
                ite = iso_Padova->interp_es(Z,tau,M,req_bands);
            else if(which=="Dartmouth")
                ite = iso_Dartmouth->interp_es(Z,tau,M,req_bands);
            else{
                s=0.;std::cerr<<"Isochrone type not valid"<<std::endl;
            }
            // std::cout<<ite[0]<<" "<<ite[1]<<" ";
            if(ite[0]<0. or ite[0]<Teffcut[0] or ite[0]>Teffcut[1]
               or ite[1]<loggcut[0] or ite[1]>loggcut[1])
                return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
            s=pow(10.,0.2*(Mag-ite[2])-2.); JK=ite[3]-ite[4];
            for(int i=5;i<5+mag_sf.size();++i)
                sf_mags.push_back(ite[i]+5*log10(s*100));
        }
        else{
            if(which=="BaSTI"){
                std::vector<int> near = iso->find_nearest(Z,tau,M);
                if(near[2]<0){
                    return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
                }
                double teff = iso->iso(near[0],near[1])->logTeff(near[2]);
                double logg = iso->iso(near[0],near[1])->logg(near[2]);
                // std::cout<<teff<<" "<<logg<<" ";
                if(teff<Teffcut[0] or teff>Teffcut[1] or logg<loggcut[0] or
                   logg>loggcut[1])
                        return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
                s = iso->iso(near[0],near[1])->distance(near[2], Mag, band);
                JK = iso->iso(near[0],near[1])->colour(near[2],{"J","K"});
                for(auto m: mag_sf)
                    sf_mags.push_back(iso->iso(near[0],near[1])->mag(near[2],m)+5*log10(s*100));
            }
            else if(which=="Padova"){
                std::vector<int> near = iso_Padova->find_nearest(Z,tau,M);
                if(near[2]<0){
                    return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
                }
                double teff = iso_Padova->iso(near[0],near[1])->logTeff(near[2]);
                double logg = iso_Padova->iso(near[0],near[1])->logg(near[2]);
                if(teff<Teffcut[0] or teff>Teffcut[1] or logg<loggcut[0] or
                   logg>loggcut[1])
                        return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
                s = iso_Padova->iso(near[0],near[1])->distance(near[2], Mag, band);
                JK = iso_Padova->iso(near[0],near[1])->colour(near[2],{"J","K"});
                for(auto m: mag_sf)
                    sf_mags.push_back(iso_Padova->iso(near[0],near[1])->mag(near[2],m)+5*log10(s*100));
            }
            else{
                s=0.;std::cerr<<"Isochrone type not valid"<<std::endl;
            }
        }
        if(extinct)
            JK+=EM->A_JK(l,b,s);
        if(fabs(b)<JKcut[1] and JK<JKcut[0])
            return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
        // -- Find polar coordinates and check if unbound
        VecDoub Pol = conv::GalacticToPolar({l,b,s},EDF->SunCoords());
        Pol[1]=0.;
        for(unsigned i=3;i<6;++i)Pol.push_back(f[i]);
        if(EDF->pot->H(Pol)>0){
            return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
        }

        // -- Calculate log-likelihood
        // -- Jacobian for d3x <-> dldbdI
        // -- s^2 cos(b) for d3x <-> dldbds and s/5. for ds <-> dI
        double Jac = s*s*s*cos(b)/5.;
        bool wf = true;
        VecDoub Actions = EDF->ActionCalculator->actions(Pol,&wf);
        double edf=log(EDF->fullDF_actions_Z(Actions, tau, Z));
        double ll = log(Jac*Mass_Function->MF(M, tau, Z))+edf;

        if(std::isnan(ll) or std::isinf(ll) or ll!=ll){
            return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
        }
        if(extinct){
            double A=0.;
            A = EM->A_V(l,b,s);
            Mag+=A*EM->extinct_const(band);
            for(int i=0;i<sf_mags.size();++i){
                sf_mags[i]+=A*EM->extinct_const(mag_sf[i]);
            }
        }
        // -- Apply selection function
        double sf = RSF->evaluate(l,b,sf_mags);
        if(sf>0) ll += log(sf);
        else{
            return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
        }
        return {ll,Actions[0],Actions[1],Actions[2],Actions[3],Actions[4],Actions[5],edf};
    }

    double LogL_sample_lb(boost::python::numeric::array f,std::string which="BaSTI", std::string band="I", bool extinct=true,bool interp=false, double modbcut=0.,double deccut=PI,VecDoub JKcut={0.,0.}, VecDoub Teffcut = {-1e6,1e6}, VecDoub loggcut = {-1e6,1e6},VecDoub FeHcut={-1e6,1e6}){
        //
        // Calculates log-likelihood given f
        // f = tau, Z, M, vR, vp,  vz, Mag, l, b
        //
        // Cut in modb in radians
        // deccut>0 we cut everything with dec<deccut-PI/2
        // deccut<0 we cut everything with dec>deccut+PI/2
        // JKcut = {Cut value (i.e. keep everything greater than cut),
        //          only apply cut for stars with fabs(b)<JKcut[1]}
        // Teffcut = Keep all Teffcut[0]<log10(Teff)<Teffcut[1]
        // loggcut = Keep all loggcut[0]<logg<loggcut[1]
        // FeHcut = Keep all FeHcut[0]<Z<FeHcut[1]
        VecDoub fvec(9,0.);
        for(unsigned i=0;i<9;++i)
            fvec[i]=extract<double>(f[i]);
        return LogL_sample(&fvec[0],which, band, extinct, interp, modbcut, deccut, JKcut, Teffcut, loggcut, FeHcut)[0];
    }
    double LogL_sample_radec_vec(double *f,std::string which="BaSTI", std::string band="I", bool extinct=true,bool interp=false, double modbcut=0.,double deccut=PI,VecDoub JKcut={0.,0.}, VecDoub Teffcut = {-1e6,1e6}, VecDoub loggcut = {-1e6,1e6},VecDoub FeHcut={-1e6,1e6}){
        //
        // Calculates log-likelihood given f
        // f = tau, Z, M, vR, vp,  vz, Mag, ra, dec
        //
        // Cut in modb in radians
        // deccut>0 we cut everything with dec<deccut-PI/2
        // deccut<0 we cut everything with dec>deccut+PI/2
        // JKcut = {Cut value (i.e. keep everything greater than cut),
        //          only apply cut for stars with fabs(b)<JKcut[1]}
        // Teffcut = Keep all Teffcut[0]<log10(Teff)<Teffcut[1]
        // loggcut = Keep all loggcut[0]<logg<loggcut[1]
        // FeHcut = Keep all FeHcut[0]<Z<FeHcut[1]
        if(f[7]<0. or f[7]>2.*PI or f[8]<-PI/2. or f[8]>PI/2.)
            return -std::numeric_limits<double>::infinity();
        VecDoub gal = conv::EquatorialToGalactic({f[7],f[8],1.});
        double Jac = log(cos(f[8]))-log(cos(gal[1]));
        f[7]=gal[0];f[8]=gal[1];
        return LogL_sample(f,which, band, extinct, interp, modbcut, deccut, JKcut, Teffcut, loggcut, FeHcut)[0]+Jac;
    }
    double LogL_sample_radec(VecDoub f,std::string which="BaSTI", std::string band="I", bool extinct=true,bool interp=false, double modbcut=0.,double deccut=PI,VecDoub JKcut={0.,0.}, VecDoub Teffcut = {-1e6,1e6}, VecDoub loggcut = {-1e6,1e6},VecDoub FeHcut={-1e6,1e6}){
        //
        // Calculates log-likelihood given f
        // f = tau, Z, M, vR, vp,  vz, Mag, ra, dec
        //
        // Cut in modb in radians
        // deccut>0 we cut everything with dec<deccut-PI/2
        // deccut<0 we cut everything with dec>deccut+PI/2
        // JKcut = {Cut value (i.e. keep everything greater than cut),
        //          only apply cut for stars with fabs(b)<JKcut[1]}
        // Teffcut = Keep all Teffcut[0]<log10(Teff)<Teffcut[1]
        // loggcut = Keep all loggcut[0]<logg<loggcut[1]
        // FeHcut = Keep all FeHcut[0]<Z<FeHcut[1]
        if(f[7]<0. or f[7]>2.*PI or f[8]<-PI/2. or f[8]>PI/2.)
            return -std::numeric_limits<double>::infinity();
        VecDoub gal = conv::EquatorialToGalactic({f[7],f[8],1.});
        double Jac = log(cos(f[8]))-log(cos(gal[1]));
        f[7]=gal[0];f[8]=gal[1];
        return LogL_sample(&f[0],which, band, extinct, interp, modbcut, deccut, JKcut, Teffcut, loggcut, FeHcut)[0]+Jac;
    }
    VecDoub LogL_sample_radec_full(VecDoub f,std::string which="BaSTI", std::string band="I", bool extinct=true,bool interp=false, double modbcut=0.,double deccut=PI,VecDoub JKcut={0.,0.}, VecDoub Teffcut = {-1e6,1e6}, VecDoub loggcut = {-1e6,1e6},VecDoub FeHcut={-1e6,1e6}){
        //
        // Calculates log-likelihood given f
        // f = tau, Z, M, vR, vp,  vz, Mag, ra, dec
        //
        // Cut in modb in radians
        // deccut>0 we cut everything with dec<deccut-PI/2
        // deccut<0 we cut everything with dec>deccut+PI/2
        // JKcut = {Cut value (i.e. keep everything greater than cut),
        //          only apply cut for stars with fabs(b)<JKcut[1]}
        // Teffcut = Keep all Teffcut[0]<log10(Teff)<Teffcut[1]
        // loggcut = Keep all loggcut[0]<logg<loggcut[1]
        // FeHcut = Keep all FeHcut[0]<Z<FeHcut[1]
        //
        // returns ll, sf, Polar coords
        if(f[7]<0. or f[7]>2.*PI or f[8]<-PI/2. or f[8]>PI/2.)
            return {-std::numeric_limits<double>::infinity(),0.,0.,0.,0.,0.,0.,0.};
        VecDoub gal = conv::EquatorialToGalactic({f[7],f[8],1.});
        double Jac = log(cos(f[8]))-log(cos(gal[1]));
        f[7]=gal[0];f[8]=gal[1];
        VecDoub result = LogL_sample(&f[0],which, band, extinct, interp, modbcut, deccut, JKcut, Teffcut, loggcut, FeHcut);
        result[0]+=Jac;
        return result;
    }
    // double LogL_sample_radec(boost::python::numeric::array f,std::string which="BaSTI", std::string band="I", bool extinct=true,bool interp=false, double modbcut=0.,double deccut=PI,VecDoub JKcut={0.,0.}, VecDoub Teffcut = {-1e6,1e6}, VecDoub loggcut = {-1e6,1e6}){
    //     //
    //     // Calculates log-likelihood given f
    //     // f = tau, Z, M, vR, vp,  vz, Mag, ra, dec
    //     //
    //     // Cut in modb in radians
    //     // deccut>0 we cut everything with dec<deccut-PI/2
    //     // deccut<0 we cut everything with dec>deccut+PI/2
    //     // JKcut = {Cut value (i.e. keep everything greater than cut),
    //     //          only apply cut for stars with fabs(b)<JKcut[1]}
    //     // Teffcut = Keep all Teffcut[0]<log10(Teff)<Teffcut[1]
    //     // loggcut = Keep all loggcut[0]<logg<loggcut[1]
    //     VecDoub fvec(9,0.);
    //     for(unsigned i=0;i<9;++i)
    //         fvec[i]=extract<double>(f[i]);
    //     return LogL_sample_radec_vec(&f[0],which, band, extinct, interp, modbcut, deccut, JKcut, Teffcut, loggcut);
    // }

    double chemDF_actions(VecDoub Acts_Z){
        return EDF->chemDF_actions(Acts_Z,Acts_Z[6]);
    }
    double log_df(VecDoub Pol, double tau, double Z){
            return log(EDF->full_DF_Z(Pol,tau,Z));
        }
    double log_dfprior(VecDoub Pol, double tau, double Z){
        return log(prior->prior_polar(Pol,Z,tau));
    }
    double log_df_actions(VecDoub Pol, VecDoub Acts, double tau, double Z){
        return log(EDF->fullDF_actions_Z(Acts,tau,Z));
    }
    double evaluate_sf(double l, double b, VecDoub sf_mags){
        return RSF->evaluate(l,b,sf_mags);
    }

    double check_sf(double tau, double Z, double M, double l, double b, VecDoub m,std::string which,bool extinct=true){

        if(l<0. or l>2.*PI or b<-PI/2. or b>PI/2.) return 0.;
        std::vector<int> near;
        if(which=="BaSTI")
            near = iso->find_nearest(Z,tau,M);
        else if(which=="Padova")
            near = iso_Padova->find_nearest(Z,tau,M);
        if(near[2]<0 and M>0.){
            return 0.;
        }
        int n=0;
        if(extinct){
            auto mbands =RSF->mag_band();
            double s;
            if(which=="BaSTI")
                s = iso->iso(near[0],near[1])->distance(near[2], m[0], mbands[0]);
            else if(which=="Padova")
                s = iso_Padova->iso(near[0],near[1])->distance(near[2], m[0], mbands[0]);
            else{
                s=0.;std::cerr<<"Isochrone type not valid"<<std::endl;
            }
            double A=0.;
            A = EM->A_V(l,b,s);
            for(int i=0;i<m.size();++i){
                m[i]+=A*EM->extinct_const(mbands[i]);
            }
        }
        return RSF->evaluate(l,b,m);
    }

    double check_color(double tau, double Z, double M, std::vector<std::string> band, VecDoub colorlimits,std::string which){
        double s, color;
        if(which=="BaSTI"){
            std::vector<int> near = iso->find_nearest(Z,tau,M);
            if(near[2]<0){
                return true;
            }
            color = iso->iso(near[0],near[1])->mag(near[2],band[0])-iso->iso(near[0],near[1])->mag(near[2],band[1]);
        }
        else if(which=="Padova"){
            std::vector<int> near = iso_Padova->find_nearest(Z,tau,M);
            if(near[2]<0){
                return false;
            }
            color = iso_Padova->iso(near[0],near[1])->mag(near[2],band[0])-iso_Padova->iso(near[0],near[1])->mag(near[2],band[1]);
        }
        if(color<colorlimits[0] or color>colorlimits[1]) return false;
        else return true;
    }

    double los_magbox_LogL_sample(boost::python::numeric::array f, double ll, double bb, std::string band, VecDoub limits, std::string which, bool RcPorZ=false, VecDoub JminusKcut={-10.,10.},VecDoub loggcut={-2000.,2000.}, VecDoub Teffcut = {.01,20000.}, VecDoub extinctflags={1.,0.}, double fieldradius=-10.,bool interp=false,VecDoub break_par={1000.,1000.},VecDoub TgJK_err={0.,0.,0.}){
        //
        // Calculates log-likelihood given f
        // f = tau, Rc/Z, M, vR, vp,  vz, apparent magnitude (V),
        //      and l and b if fieldradius>0
        bool extinct = extinctflags[0]!=0.;
        bool dered = extinctflags[1]!=0.;

        double lower_limit=limits[0];
        double upper_limit=limits[1];

    	double l=0.,b=0.;
        if(fieldradius>0.){
            if(len(f)>7){
                l = extract<double>(f[7]);
                b = extract<double>(f[8]);
                if(pow(cos(bb)*(l-ll),2.)+pow(b-bb,2.)>pow(fieldradius,2.))
                    return -std::numeric_limits<double>::infinity();
            }
            else std::cerr<<"For sampling l and b must specify fieldradius (in radians) and ll and bb are field centre\n";
        }
        else{l=ll;b=bb;}
        // -- Extract quantities from f and check if allowed
        double tau = extract<double>(f[0]);
        if(tau>=EDF->get_taum() or tau<0.)
            return -std::numeric_limits<double>::infinity();

        double ZorRcP = extract<double>(f[1]), Z, RcP;
        if(RcPorZ==false){
            Z = ZorRcP;
            if(Z<-3. or Z>1.) return -std::numeric_limits<double>::infinity();
        }
        if(RcPorZ==true){
            RcP=RcPorZ;
            if(RcP<0.) return -std::numeric_limits<double>::infinity();
            Z = EDF->FeH(tau, RcP);
        }
        double M = extract<double>(f[2]);
        if(M<0.) return -std::numeric_limits<double>::infinity();

        double V = extract<double>(f[6]);
        if(V>upper_limit and !extinct)
            return -std::numeric_limits<double>::infinity();
        if(V<lower_limit and !extinct)
            return -std::numeric_limits<double>::infinity();

        // -- Find nearest isochrone point or interp
        double s=0.,JK=0.;
        if(interp){
            VecDoub ite;
            if(which=="BaSTI")
                ite = iso->interp_es(Z,tau,M,{band,"J","K"});
            else if(which=="Padova")
                ite = iso_Padova->interp_es(Z,tau,M,{band,"J","K"});
            else if(which=="Dartmouth")
                ite = iso_Dartmouth->interp_es(Z,tau,M,{band,"J","K"});
            else{
                s=0.;std::cerr<<"Isochrone type not valid"<<std::endl;
            }
            ite[0]+=rnGauss->nextnumber()*TgJK_err[0];
            ite[1]+=rnGauss->nextnumber()*TgJK_err[1];
            if(ite[0]<log10(Teffcut[0]) or ite[0]>log10(Teffcut[1]))
                return -std::numeric_limits<double>::infinity();
            if(ite[1]<loggcut[0] or ite[1]>loggcut[1])
                return -std::numeric_limits<double>::infinity();
            JK=ite[3]-ite[4];JK+=rnGauss->nextnumber()*TgJK_err[2];
            s=pow(10.,0.2*(V-ite[2])-2.);
            if(extinct and !dered)
                JK+=EM->A_JK(l,b,s);
            if(JK<JminusKcut[0] or JK>JminusKcut[1])
                return -std::numeric_limits<double>::infinity();
        }
        else{
            std::vector<int> near; double lg=0., Teff=0.;
            if(which=="BaSTI"){
                near = iso->find_nearest(Z,tau,M);
                lg = iso->iso(near[0],near[1])->logg(near[2]);
                Teff = iso->iso(near[0],near[1])->logTeff(near[2]);
                s = iso->iso(near[0],near[1])->distance(near[2], V, band);
                JK = iso->iso(near[0],near[1])->colour(near[2],{"J","K"});
            }
            else if(which=="Padova"){
                near = iso_Padova->find_nearest(Z,tau,M);
                lg = iso_Padova->iso(near[0],near[1])->logg(near[2]);
                Teff = iso_Padova->iso(near[0],near[1])->logTeff(near[2]);
                s = iso_Padova->iso(near[0],near[1])->distance(near[2], V, band);
                JK = iso_Padova->iso(near[0],near[1])->colour(near[2],{"J","K"});
            }
            else if(which=="Dartmouth"){
                near = iso_Dartmouth->find_nearest(Z,tau,M);
                lg = iso_Dartmouth->iso(near[0],near[1])->logg(near[2]);
                Teff = iso_Dartmouth->iso(near[0],near[1])->logTeff(near[2]);
                s = iso_Dartmouth->iso(near[0],near[1])->distance(near[2], V, band);
                JK = iso_Dartmouth->iso(near[0],near[1])->colour(near[2],{"J","K"});
            }
            else{
                s=0.;std::cerr<<"Isochrone type not valid"<<std::endl;
            }
            if(near[2]<0)
                return -std::numeric_limits<double>::infinity();

            Teff+=rnGauss->nextnumber()*TgJK_err[0];
            lg+=rnGauss->nextnumber()*TgJK_err[1];
            JK+=rnGauss->nextnumber()*TgJK_err[2];
            if(lg<loggcut[0] or lg>loggcut[1])
                return -std::numeric_limits<double>::infinity();
            if(Teff<log10(Teffcut[0]) or Teff>log10(Teffcut[1]))
                return -std::numeric_limits<double>::infinity();
            if(extinct and !dered)
                JK+=EM->A_JK(l,b,s);
            if(JK<JminusKcut[0] or JK>JminusKcut[1])
                return -std::numeric_limits<double>::infinity();
        }
        if(extinct){
            double A=0.;
            if(band=="V") A = EM->A_V(l,b,s);
            if(band=="H") A = EM->A_H(l,b,s);
            if(band=="J") A = EM->A_J(l,b,s);
            V+=A;
        }
        if(V>upper_limit and extinct)
            return -std::numeric_limits<double>::infinity();
        if(V<lower_limit and extinct)
            return -std::numeric_limits<double>::infinity();
	    // -- Find polar coordinates and check if unbound
        VecDoub Pol = conv::GalacticToPolar({l,b,s},EDF->SunCoords());
        for(unsigned i=3;i<6;++i)Pol.push_back(extract<double>(f[i]));
        Pol[1]=0.;
        if(EDF->pot->H(Pol)>0.)
            return -std::numeric_limits<double>::infinity();
        // -- Calculate log-likelihood
        // -- Jacobian for d3x <-> dldbdI
        // -- s^2 cos(b) for d3x <-> dldbds and s/5./log10e for ds <-> dI
        double Jac = s*s*s*cos(b)/11.51;
        double loglike;
        if(RcPorZ==true)
            loglike = log(Jac*EDF->full_DF(Pol, tau, RcP)*Mass_Function->MF(M, tau, Z));
        else
            loglike = log(Jac*EDF->full_DF_Z(Pol, tau, Z)*Mass_Function->MF(M, tau, Z));
        if(std::isnan(loglike) or std::isinf(loglike) or loglike!=loglike)
            return -std::numeric_limits<double>::infinity();
        if(V>break_par[0]+1./break_par[1])
    		return -std::numeric_limits<double>::infinity();
        if(V>break_par[0])
    		loglike+=log(1.-break_par[1]*(V-break_par[0]));
        return loglike;
    }

    VecDoub get_extra_data(boost::python::numeric::array f, std::string band, std::string which, bool RcPorZ=false,bool extinct=true,bool interp=false){
        // Take a sample f and calculate the additional variables from the
        // model
        double tau = extract<double>(f[0]);
        double ZorRcP = extract<double>(f[1]),Z, RcP;
        if(RcPorZ==false){
            Z = ZorRcP;
            RcP = EDF->RadiusFromMetal(tau,ZorRcP);
        }
        else{
            RcP = ZorRcP;
            Z = EDF->FeH(tau, ZorRcP);
        }
        double M = extract<double>(f[2]);
        double Mg = extract<double>(f[6]);
        double l = extract<double>(f[7]);
        double b = extract<double>(f[8]);

        // -- Find nearest isochrone point
        double dm=0.,s=0., JK=0., Teff=0., logg=0.,J=0.,K=0.;
        if(which=="BaSTI"){
            if(interp){
                VecDoub V = iso->interp_es(Z,tau,M,{band,"J","K"});
                Teff=V[0];logg=V[1];
                JK=V[3]-V[4];
                s=pow(10.,0.2*(Mg-V[2])-2.);
                dm = 5.*log10(100.*s);J=V[3]+dm;K=V[4]+dm;
            }
            else{
                std::vector<int> near = iso->find_nearest(Z,tau,M);
                Teff = iso->iso(near[0],near[1])->logTeff(near[2]);
                logg = iso->iso(near[0],near[1])->logg(near[2]);
                s = iso->iso(near[0],near[1])->distance(near[2], Mg, band);
                dm = 5.*log10(100.*s);
		        JK = iso->iso(near[0],near[1])->colour(near[2],{"J","K"});
                J = iso->iso(near[0],near[1])->mag(near[2],"J")+dm;
                K = iso->iso(near[0],near[1])->mag(near[2],"K")+dm;
            }
        }
        else if(which=="Padova"){
            if(interp){
                VecDoub V = iso_Padova->interp_es(Z,tau,M,{band,"J","K"});
                Teff=V[0];logg=V[1];
                JK=V[3]-V[4];
                s=pow(10.,0.2*(Mg-V[2])-2.);
                dm = 5.*log10(100.*s);J=V[3]+dm;K=V[4]+dm;
            }
            else{
                std::vector<int> near = iso_Padova->find_nearest(Z,tau,M);
                Teff = iso_Padova->iso(near[0],near[1])->logTeff(near[2]);
                logg = iso_Padova->iso(near[0],near[1])->logg(near[2]);
                s = iso_Padova->iso(near[0],near[1])->distance(near[2], Mg, band);
                dm = 5.*log10(100.*s);
        	JK = iso_Padova->iso(near[0],near[1])->colour(near[2],{"J","K"});
                J = iso_Padova->iso(near[0],near[1])->mag(near[2],"J")+dm;
                K = iso_Padova->iso(near[0],near[1])->mag(near[2],"K")+dm;
            }
        }
        else if(which=="Dartmouth"){
            if(interp){
                VecDoub V = iso_Dartmouth->interp_es(Z,tau,M,{band,"J","K"});
                Teff=V[0];logg=V[1];
                JK=V[3]-V[4];
                s=pow(10.,0.2*(Mg-V[2])-2.);
                dm = 5.*log10(100.*s);J=V[3]+dm;K=V[4]+dm;
            }
            else{
                std::vector<int> near = iso_Dartmouth->find_nearest(Z,tau,M);
                Teff = iso_Dartmouth->iso(near[0],near[1])->logTeff(near[2]);
                logg = iso_Dartmouth->iso(near[0],near[1])->logg(near[2]);
                s = iso_Dartmouth->iso(near[0],near[1])->distance(near[2], Mg, band);
                dm = 5.*log10(100.*s);
        	JK = iso_Dartmouth->iso(near[0],near[1])->colour(near[2],{"J","K"});
                J = iso_Dartmouth->iso(near[0],near[1])->mag(near[2],"J")+dm;
                K = iso_Dartmouth->iso(near[0],near[1])->mag(near[2],"K")+dm;
            }
        }
        if(extinct){
            double A=0.;
            A = EM->A_V(l,b,s)*EM->extinct_const(band);
            Mg+=A;
            if(A<0.)std::cerr<<"A<0 at l="<<l<<", b="<<b<<", s="<<s<<std::endl;
            JK+=EM->A_JK(l,b,s);
            J+=EM->A_J(l,b,s);
            K+=EM->A_K(l,b,s);
        }
        // -- Find Galactic coordinates
        VecDoub Pol = conv::GalacticToPolar({l,b,s},EDF->SunCoords());
        for(int i=3;i<6;i++)Pol.push_back(extract<double>(f[i]));
        Pol[4]*=-1.;
        VecDoub Gal = conv::PolarToGalactic(Pol,EDF->SunCoords());
        VecDoub Eq = conv::GalacticToEquatorial(Gal);
        VecDoub extras = {s,RcPorZ?Z:RcP,Teff,logg,Pol[0],Pol[1],Pol[2],
                 Gal[3],Gal[4],Gal[5],
                 Eq[0],Eq[1],Eq[4],Eq[5],JK,J,K,Mg};
        return extras;
    }

    double get_magnitude(double Z, double age, double M, std::string mag, std::string which){
        if(which=="BaSTI"){
            std::vector<int> near = iso->find_nearest(Z,age,M);
	    return iso->iso(near[0],near[1])->mag(near[2],mag);
	}
        else if(which=="Padova"){
            std::vector<int> near = iso_Padova->find_nearest(Z,age,M);
	    return iso_Padova->iso(near[0],near[1])->mag(near[2],mag);
	}
        else if(which=="Dartmouth"){
            std::vector<int> near = iso_Dartmouth->find_nearest(Z,age,M);
	    return iso_Dartmouth->iso(near[0],near[1])->mag(near[2],mag);
	}
	else{
		std::cerr<<"Invalid isochrone"<<std::endl;
		return -9999.;
	}
}

    VecDoub get_APASS(boost::python::numeric::array f, std::string band, std::string which, bool RcPorZ=false, bool extinct=true,bool interp=false){
        // Take a sample f and calculate the additional variables from the
        // model
        double tau = extract<double>(f[0]);
        double ZorRcP = extract<double>(f[1]),Z, RcP;
        if(RcPorZ==false){
            Z = ZorRcP;
            RcP = EDF->RadiusFromMetal(tau,ZorRcP);
        }
        else{
            RcP = ZorRcP;
            Z = EDF->FeH(tau, ZorRcP);
        }
        double M = extract<double>(f[2]);
        double Mg = extract<double>(f[6]);
        // -- Find nearest isochrone point
        double B=0.,V=0.,s=0.,dm=0.;
        if(which=="BaSTI"){
            if(interp){
                VecDoub ii = iso->interp_es(Z,tau,M,{band,"B","V"});
                s=pow(10.,0.2*(Mg-ii[2])-2.);
                dm = 5.*log10(100.*s);
                B=ii[3]+dm; V=ii[4]+dm;
            }
            else{
                std::vector<int> near = iso->find_nearest(Z,tau,M);
                s = iso->iso(near[0],near[1])->distance(near[2], Mg, band);
                dm = 5.*log10(100.*s);
                B = iso->iso(near[0],near[1])->mag(near[2],"B")+dm;
                V = iso->iso(near[0],near[1])->mag(near[2],"V")+dm;
            }
      	}
        else if(which=="Padova"){
            if(interp){
                VecDoub ii = iso_Padova->interp_es(Z,tau,M,{band,"B","V"});
                s=pow(10.,0.2*(Mg-ii[2])-2.);
                dm = 5.*log10(100.*s);
                B=ii[3]+dm; V=ii[4]+dm;
            }
            else{
                std::vector<int> near = iso_Padova->find_nearest(Z,tau,M);
                s = iso_Padova->iso(near[0],near[1])->distance(near[2], Mg, band);
                dm = 5.*log10(100.*s);
                B = iso_Padova->iso(near[0],near[1])->mag(near[2],"B")+dm;
                V = iso_Padova->iso(near[0],near[1])->mag(near[2],"V")+dm;
            }
        }
      	else std::cerr<<"Isochrone system not valid"<<std::endl;
        if(extinct){
      	    double l = extract<double>(f[7]);
      	    double b = extract<double>(f[8]);
                  B += EM->A_B(l,b,s);
                  V += EM->A_V(l,b,s);
      	}
      	return {B,V};
    }


    double get_EBV(double l, double b, double s){
    	return EM->A_V(l,b,s)/3.1;
    }

    boost::python::numeric::array get_actions(boost::python::numeric::array f){
        VecDoub Pol {extract<double>(f[0])};
        for(int i=1;i<3;i++)Pol.push_back(extract<double>(f[i]));
        for(int i=3;i<6;i++)Pol.push_back(extract<double>(f[i]));
        Pol[1]=0.;
        if(EDF->pot->H(Pol)>0. or Pol[0]!=Pol[0] or std::isinf(Pol[0]) or std::isnan(Pol[0]))
            return boost::python::numeric::array(
                boost::python::make_tuple(
                 std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity()));
        bool wf = true;
        VecDoub Acts = EDF->ActionCalculator->actions(Pol,&wf);
        return boost::python::numeric::array(
                boost::python::make_tuple(
                 Acts[0],Acts[1],Acts[2],Acts[3]));
    }
    boost::python::numeric::array get_actions_freqs(boost::python::numeric::array f){
        VecDoub Pol {extract<double>(f[0])};
        for(int i=1;i<3;i++)Pol.push_back(extract<double>(f[i]));
        for(int i=3;i<6;i++)Pol.push_back(extract<double>(f[i]));
        Pol[1]=0.;
        if(EDF->pot->H(Pol)>0. or Pol[0]!=Pol[0] or std::isinf(Pol[0]) or std::isnan(Pol[0]))
            return boost::python::numeric::array(
                boost::python::make_tuple(
                 std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity(),
                 std::numeric_limits<double>::infinity()));
        bool wf = true;
        VecDoub Acts = EDF->ActionCalculator->actions(Pol,&wf);
        return boost::python::numeric::array(
                boost::python::make_tuple(
                 Acts[0],Acts[1],Acts[2],Acts[3],Acts[4],Acts[5]));
    }

    double check_highmass(double age, double Rc, std::string which="BaSTI"){
        // Checking if mass is physical at age, Z
        double Z = EDF->FeH(age, Rc);
        if(which=="BaSTI"){
            std::vector<int> near = iso->find_nearest(Z,age,0.6);
            return iso->iso(near[0],near[1])->maxmass();
        }
        else if(which=="Padova"){
            std::vector<int> near = iso_Padova->find_nearest(Z,age,0.6);
            return iso_Padova->iso(near[0],near[1])->maxmass();
        }
        else if (which=="Dartmouth"){
            std::vector<int> near = iso_Dartmouth->find_nearest(Z,age,0.6);
            return iso_Dartmouth->iso(near[0],near[1])->maxmass();
        }
        else{std::cerr<<"Not valid isochrone grid in highmass\n";return -1.;}
    }

    double check_highmass_Z(double age, double Z, std::string which="BaSTI"){
        // Checking if mass is physical at age, Z
        if(age>EDF->MaxAge(Z) and age<12.) return -100.;
        if(which=="BaSTI"){
            std::vector<int> near = iso->find_nearest(Z,age,0.6);
            return iso->iso(near[0],near[1])->maxmass();
        }
        else if(which=="Padova"){
            std::vector<int> near = iso_Padova->find_nearest(Z,age,0.6);
            return iso_Padova->iso(near[0],near[1])->maxmass();
        }
        else if (which=="Dartmouth"){
            std::vector<int> near = iso_Dartmouth->find_nearest(Z,age,0.6);
            return iso_Dartmouth->iso(near[0],near[1])->maxmass();
        }
        else{std::cerr<<"Not valid isochrone grid in highmassZ\n";return -1.;}
    }


    double check_lowmass(double age, double Rc, std::string which="BaSTI"){
        // Checking if mass is physical at age, Rc
        double Z = EDF->FeH(age, Rc);
        if(which=="BaSTI"){
            std::vector<int> near = iso->find_nearest(Z,age,0.6);
            return iso->iso(near[0],near[1])->minmass();
        }
        else if(which=="Padova"){
            std::vector<int> near = iso_Padova->find_nearest(Z,age,0.6);
            return iso_Padova->iso(near[0],near[1])->minmass();
        }
        else if (which=="Dartmouth"){
            std::vector<int> near = iso_Dartmouth->find_nearest(Z,age,0.6);
            return iso_Dartmouth->iso(near[0],near[1])->minmass();
        }
        else{std::cerr<<"Not valid isochrone grid in lowmass\n";return -1.;}
    }

    double check_lowmass_Z(double age, double Z, std::string which="BaSTI"){
        // Checking if mass is physical at age, Z
        if(age>EDF->MaxAge(Z) and age<12.) return -100.;
        if(which=="BaSTI"){
            std::vector<int> near = iso->find_nearest(Z,age,0.6);
            return iso->iso(near[0],near[1])->minmass();
        }
        else if(which=="Padova"){
            std::vector<int> near = iso_Padova->find_nearest(Z,age,0.6);
            return iso_Padova->iso(near[0],near[1])->minmass();
        }
        else if (which=="Dartmouth"){
            std::vector<int> near = iso_Dartmouth->find_nearest(Z,age,0.6);
            return iso_Dartmouth->iso(near[0],near[1])->minmass();
        }
        else{std::cerr<<"Not valid isochrone grid in lowmassZ\n";return -1.;}
    }

    double minZ(void){
        return EDF->minZ();
    }

    double maxZ(void){
        return EDF->maxZ();
    }

    bool check_radius_positive(double age, double Z){
        if((age>EDF->MaxAge(Z) or EDF->RadiusFromMetal(age,Z)<0.)and age<11.) return false;
        else return true;
    }

    bool check_JK_logg_cut(double age, double Z, double M, double b, double bcut, VecDoub JminusKcut,VecDoub loggcut,VecDoub Teffcut, std::string which="BaSTI", bool interp=false){
        double JK,lg,Teff;
        if(which=="BaSTI"){
            if(interp){
                VecDoub ite = iso->interp_es(Z,age,M,{"J","K"});
                Teff=ite[0];lg=ite[1];JK=ite[2]-ite[3];
            }
            else{
                std::vector<int> near = iso->find_nearest(Z,age,M);
                JK = iso->iso(near[0],near[1])->colour(near[2],{"J","K"});
                lg = iso->iso(near[0],near[1])->logg(near[2]);
                Teff = iso->iso(near[0],near[1])->logTeff(near[2]);
            }
        }
	else if(which=="Padova"){
            if(interp){
                VecDoub ite = iso_Padova->interp_es(Z,age,M,{"J","K"});
                Teff=ite[0];lg=ite[1];JK=ite[2]-ite[3];
            }
            else{
                std::vector<int> near = iso_Padova->find_nearest(Z,age,M);
                JK = iso_Padova->iso(near[0],near[1])->colour(near[2],{"J","K"});
                lg = iso_Padova->iso(near[0],near[1])->logg(near[2]);
                Teff = iso_Padova->iso(near[0],near[1])->logTeff(near[2]);
            }
        }
        else if (which=="Dartmouth"){
            if(interp){
                VecDoub ite = iso_Dartmouth->interp_es(Z,age,M,{"J","K"});
                Teff=ite[0];lg=ite[1];JK=ite[2]-ite[3];
            }
            else{
                std::vector<int> near = iso_Dartmouth->find_nearest(Z,age,M);
                JK = iso_Dartmouth->iso(near[0],near[1])->colour(near[2],{"J","K"});
                lg = iso_Dartmouth->iso(near[0],near[1])->logg(near[2]);
                Teff = iso_Dartmouth->iso(near[0],near[1])->logTeff(near[2]);
            }
        }
        else{
            std::cerr<<"Invalid isochrone type in check_JK_cut\n";
            return false;
        }
        if((fabs(b)<=bcut*PI/180. and (JK<JminusKcut[0] or JK>JminusKcut[1])) or lg>loggcut[1] or lg<loggcut[0] or Teff<log10(Teffcut[0]) or Teff>log10(Teffcut[1]))
            return false;
        else return true;
    }

    VecDoub process_data(VecDoub Eq){
        // printVector(Eq);
        VecDoub Gal=conv::EquatorialToGalactic(Eq);
        // printVector(Gal);
        VecDoub Pol=conv::GalacticToPolar(Gal,EDF->SunCoords());
        VecDoub Acts(4,0);
        double phi = Pol[1];
        Pol[4]*=-1.;Pol[1]=0.;
        //if(EDF->pot->H(Pol)>0. or fabs(Pol[4])>500.)
        if(EDF->pot->H(Pol)>0.)
            for(int i=0;i<4;++i)Acts[i]=std::numeric_limits<double>::infinity();
        else{
            bool wf = true;
            Acts = EDF->ActionCalculator->actions(Pol,&wf);
        }
        return {Gal[0],Gal[1],Gal[2],Gal[3],Gal[4],Gal[5],
            Pol[0],phi,Pol[2],Pol[3],Pol[4],Pol[5],
            Acts[0],Acts[1],Acts[2],Acts[3]};
    }

    double check_dec(double ll, double bb){
        // returns dec
        VecDoub eq = conv::GalacticToEquatorial({ll,bb,1.});
        return eq[1];
    }
    double get_extinct(double l, double b, double s, std::string band){
        return EM->A_V(l,b,s)*EM->extinct_const(band);
    }
    double get_extinct_H(double l, double b, double s){
        return EM->A_H(l,b,s);
    }
    double get_extinct_J(double l, double b, double s){
        return EM->A_J(l,b,s);
    }
    double get_extinct_I(double l, double b, double s){
        return EM->A_I(l,b,s);
    }
    void turn_on_halo(void){
        EDF->TurnOnHalo();
    }
    void turn_off_halo(void){
        EDF->TurnOffHalo();
    }
    void print_SF(VecDoub mag){
        for(double l=0.0001;l<2.*PI;l+=2.*PI/100.)
        for(double b=-PI/2.+0.0001;b<PI/2.;b+=PI/50.)
            std::cout<<l<<" "<<b<<" "<<RSF->evaluate(l,b,mag)<<std::endl;
    }

    struct norm_st: edf_norm_struct{
        std::string which="BaSTI";
        std::string band = "I";
        bool extinct=true;
        bool interp=false;
        double modbcut=0.;
        double deccut=PI;
        VecDoub JKcut={0.,0.};
        VecDoub Teffcut={-1e6,1e6};
        VecDoub loggcut={-1e6,1e6};
        norm_st(VecDoub x2min,
                VecDoub x2max,
                std::string which="BaSTI",
                std::string band = "I",
                  bool extinct=true,
                  bool interp=false,
                  double modbcut=0.,
                  double deccut=PI,
                  VecDoub JKcut={0.,0.},
                  VecDoub Teffcut={-1e6,1e6},
                  VecDoub loggcut={-1e6,1e6})
            :edf_norm_struct(nullptr,x2min,x2max,0.), which(which), band(band),
                extinct(extinct), interp(interp), modbcut(modbcut),
                deccut(deccut),JKcut(JKcut),Teffcut(Teffcut),loggcut(loggcut){};
    };

    static int norm_integrand_cuba(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
        norm_st *P = (norm_st *) fdata;
        VecDoub y2(9,0.);
        for(unsigned int i=0;i<P->x2min.size();i++){
            y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
            std::cout<<(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i]<<" ";
        }
        double rslt=LogL_sample_radec(y2,P->which,P->band,P->extinct,
                                          P->interp, P->modbcut,P->deccut,
                                          P->JKcut,P->Teffcut, P->loggcut);
        if(rslt!=rslt)
            fval[0]=0.;
        else
            fval[0] = exp(rslt);
        std::cout<<fval[0]<<std::endl;
        return 0;
    }

    double find_normalization(std::string which="BaSTI",
                              bool extinct=true,
                              bool interp=false,
                              double modbcut=0.,
                              double deccut=PI,
                              VecDoub JKcut={0.,0.},
                              VecDoub Teffcut={-1e6,1e6},
                              VecDoub loggcut={-1e6,1e6}){
        // tau, Z, M, vR, vp,  vz, I, ra, dec
        double tmin=0.;
        double tmax=EDF->get_taum();
        double Zmax=EDF->get_max_Z();
        double Zmin=EDF->get_min_Z();
        double minMass = 0.8, maxMass = 4.;
        double minMag = RSF->min_mag();
        double maxMag = RSF->max_mag();
        double ramin = 0.;
        double ramax = 2.*PI;
        double decmin = -.5*PI;
        double decmax = .5*PI;

        double vmin = -300.;
        double vmax =  300.;
        double vphimin = -300.;
        double vphimax =  400.;

        if(deccut<PI){
            if(deccut>0.)
                decmin = deccut-.5*PI;
            else if(deccut<0.)
                decmax = deccut+.5*PI;
        }

        std::vector<VecDoub> peaks = {{1.,0.,1.,10.,0.,220.,0.,1.5*PI,-.25*PI}};

        minMag = 10.;
        maxMag = 10.2;
        // ramin = 1.5*PI;
        // ramax = 1.55*PI;
        // decmin = -.25*PI;
        // decmax = -.24*PI;

        // vmin = -30.;
        // vmax = -20.;
        // vphimin = 200.;
        // vphimax = 210.;

        VecDoub x2min = {tmin,Zmin,minMass,
                         vmin, vphimin, vmin,
                         minMag, ramin, decmin};
        VecDoub x2max = {tmax,Zmax,maxMass,
                          vmax, vphimax, vmax,
                          maxMag, ramax, decmax};
        double IE = 1e-3;
        double err;
        norm_st P(x2min,x2max,which,RSF->mag_band()[0],extinct,interp,modbcut,
                  deccut,JKcut,Teffcut,loggcut);
        // factor of two for velocity space
        return 2.*edf_integrate(norm_integrand_cuba,&P,IE,0,"Divonne",&err,0,&peaks);
    }
    VecDoub get_limits(void){
        double tmin=0.;
        double tmax=EDF->get_taum();
        double Zmax=EDF->get_max_Z();
        double Zmin=EDF->get_min_Z();
        double minMass = 0.8, maxMass = 4.;
        double minMag = RSF->min_mag();
        double maxMag = RSF->max_mag();
        double lmin = 0.;
        double lmax = 2.*PI;
        double bmin = -.5*PI;
        double bmax = .5*PI;

        double vmin = -300.;
        double vmax =  300.;
        double vphimin = -300.;
        double vphimax =  400.;
        return {tmin,tmax,Zmin,Zmax,minMass,maxMass,vmin,vmax,vphimin,vphimax,
                vmin,vmax,minMag,maxMag,lmin,lmax,bmin,bmax};
    }

    void load_mc_samples(std::string fname){
        H5File fin(fname,H5F_ACC_RDONLY);
        hdf5_read_1D_vector(fin,mc_samples,"mc_samples");
        fin.close();
    }
    double eval_mc_integral(int nstart=0, int nend=1){
        int ndim=9; VecDoub J(6,0.);
        double LL=0., tau, Z, ledf, edf;
        for(int i=nstart;i<nend;++i){
            for(unsigned k=0;k<6;++k)
                J[k]=mc_samples[i*ndim+2+k];
            tau =mc_samples[i*ndim];
            Z = mc_samples[i*ndim+1];
            ledf=mc_samples[i*ndim+ndim-1];
            edf=EDF->fullDF_actions_Z(J,tau,Z);
            if(edf>0.) LL+=exp(log(edf)-ledf);
        }
        return LL;
    }

    int load_data_samples(std::string fname, bool set_med=false){
        H5File fin(fname,H5F_ACC_RDONLY);
        hdf5_read_1D_vector(fin,data_samples,"data_samples");
        fin.close();
        if(set_med){
            int ndim=14;
            VecDoub med(12,0.);
            unsigned datalen=data_samples.size()/ndim;
            for(int i=0;i<datalen;++i){
                for(int k=0;k<6;++k)
                    med[k]+=data_samples[i*ndim+k];
                for(int k=8;k<ndim;++k)
                    med[k-2]+=data_samples[i*ndim+k];
            }
            for(int k=0;k<12;++k)
                med[k]/=(double)(datalen);
            for(int i=0;i<datalen;++i){
                for(int k=0;k<6;++k)
                    data_samples[i*ndim+k]=med[k];
                for(int k=8;k<ndim;++k)
                    data_samples[i*ndim+k]=med[k-2];
            }
        }
        return data_samples.size();
    }
    double logsumexp(double *nums, size_t ct) {
      double max_exp = nums[0], sum = 0.0;
      size_t i;
      for (i = 1 ; i < ct ; i++)
        if (nums[i] > max_exp)
          max_exp = nums[i];
      for (i = 0; i < ct ; i++)
        sum += exp(nums[i] - max_exp);
      return log(sum) + max_exp;
    }
    double eval_data_likelihood(int nstart=0, int nend=1, int Nsamples=100, bool prior_weights=true){
        double lNsamples=log(Nsamples);
        unsigned ndim=15;
        unsigned nstride=ndim*Nsamples;
        unsigned skip=0;
        VecDoub err_s(Nsamples,0.), J(6,0.), pol(6,0.);
        double logL=0.,tau,Z,prior;
        for(unsigned i=nstart;i<nend;++i){
            for(unsigned j=0;j<Nsamples;++j){
                skip=i*nstride+j*ndim;
                for(unsigned k=0;k<6;++k){
                    pol[k]=data_samples[skip+k];
                    J[k]=data_samples[skip+8+k];
                }
                tau=data_samples[skip+6];Z=data_samples[skip+7];
                // Prior weights include the Jacobian factor
                //prior=data_samples[skip+14];
                err_s[j]=log_df_actions(pol,J,tau,Z);
    	        if(std::isinf(err_s[j]) or std::isnan(err_s[j]) or err_s[j]!=err_s[j])
                    err_s[j]=-1000.;
                if(prior_weights)
		    err_s[j]-=data_samples[skip+14];
            }
            logL+=logsumexp(&err_s[0],Nsamples)-lNsamples;
	}
        return logL;
    }
}


BOOST_PYTHON_MODULE(edf_sampling)
{
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
    def("get_params", edf_sampling::get_params);
    def("set_params", edf_sampling::set_params);
    def("set_params_json", edf_sampling::set_params_json);
    def("get_Vc", edf_sampling::get_Vc);
    def("LogL_sample", edf_sampling::LogL_sample_lb);
    def("LogL_sample_radec", edf_sampling::LogL_sample_radec);
    def("LogL_sample_radec_full", edf_sampling::LogL_sample_radec_full);
    def("los_magbox_LogL_sample", edf_sampling::los_magbox_LogL_sample);
    def("setup", edf_sampling::setup);
    def("load_selection_function", edf_sampling::load_selection_function);
    def("load_selection_function_stack",
        edf_sampling::load_selection_function_stack);
    def("init_isochrone", edf_sampling::init_isochrone);
    def("print_SF", edf_sampling::print_SF);
    def("load_emap", edf_sampling::load_emap);
    def("minZ", edf_sampling::minZ);
    def("maxZ", edf_sampling::maxZ);
    def("chemDF_actions", edf_sampling::chemDF_actions);
    def("get_extra_data", edf_sampling::get_extra_data);
    def("get_APASS", edf_sampling::get_APASS);
    def("get_magnitude", edf_sampling::get_magnitude);
    def("get_actions", edf_sampling::get_actions);
    def("get_actions_freqs", edf_sampling::get_actions_freqs);
    def("check_color", edf_sampling::check_color);
    def("check_highmass",edf_sampling::check_highmass);
    def("check_highmass_Z",edf_sampling::check_highmass_Z);
    def("check_lowmass",edf_sampling::check_lowmass);
    def("check_lowmass_Z",edf_sampling::check_lowmass_Z);
    def("check_radius_positive",edf_sampling::check_radius_positive);
    def("check_JK_logg_cut",edf_sampling::check_JK_logg_cut);
    def("check_dec",edf_sampling::check_dec);
    def("check_sf",edf_sampling::check_sf);
    def("get_extinct",edf_sampling::get_extinct);
    def("get_extinct_H",edf_sampling::get_extinct_H);
    def("get_extinct_J",edf_sampling::get_extinct_J);
    def("get_extinct_I",edf_sampling::get_extinct_I);
    def("get_EBV",edf_sampling::get_EBV);
    def("process_data",edf_sampling::process_data);
    def("turn_off_halo",edf_sampling::turn_off_halo);
    def("turn_on_halo",edf_sampling::turn_on_halo);
    def("find_normalization",edf_sampling::find_normalization);
    def("get_limits",edf_sampling::get_limits);
    def("log_df",edf_sampling::log_df);
    def("log_dfprior",edf_sampling::log_dfprior);
    def("log_df_actions",edf_sampling::log_df_actions);
    def("evaluate_sf",edf_sampling::evaluate_sf);
    def("load_mc_samples",edf_sampling::load_mc_samples);
    def("eval_mc_integral",edf_sampling::eval_mc_integral);
    def("load_data_samples",edf_sampling::load_data_samples);
    def("eval_data_likelihood",edf_sampling::eval_data_likelihood);
    class_<isochrone>("isochrone",no_init)
    .def("initial_mass",&isochrone::initial_mass)
    .def("mass",&isochrone::mass)
    .def("logg",&isochrone::logg)
    .def("logTeff",&isochrone::logTeff);
    class_<isochrone_johnson, bases<isochrone> >("isochrone_johnson",no_init)
    .def("fill",&isochrone_johnson::fill);
      to_python_converter<VecDoub, vector_to_ndarray<double>>();
      vector_from_ndarray<double>();
      to_python_converter<std::vector<std::string>, vector_to_ndarray<std::string>>();
      vector_from_ndarray<std::string>();
      import_array();
}
