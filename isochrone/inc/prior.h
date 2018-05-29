#include "utils_iso.h"

class galaxy_prior{
    protected:
        bool eval_prior;
        VecDoub SolarPosition;
        double R0, z0;
    public:
        galaxy_prior(VecDoub SolarPosition, bool ep=false):eval_prior(ep),SolarPosition(SolarPosition){R0=SolarPosition[0];z0=SolarPosition[1];}
        VecDoub solar(void){return SolarPosition;}
        virtual double prior(VecDoub x, double Z, double tau){
            return 1.;
        }
        virtual double prior(VecDoub x, double Z, double tau, double alpha){
            return prior(x,Z,tau);
        }
        virtual double prior_polar(VecDoub Rpz, double Z, double tau){
            return prior({Rpz[0]*cos(Rpz[1]),Rpz[0]*sin(Rpz[1]),Rpz[2]},Z,tau);
        }
        virtual double prior_polar(VecDoub Rpz, double Z, double tau, double alpha){
            return prior_polar(Rpz,Z,tau);
        }
        bool bprior(void){return eval_prior;}
        void turn_on_prior(void){eval_prior=true;}
        void turn_off_prior(void){eval_prior=false;}
        void set_prior(bool bprior){eval_prior=bprior;}
};

inline double double_exponential(double R, double z, double Rd, double zd){
    // Remember need to pass abs(z) but we also normalise to solar values at z0 so pass
    // abs(z)-z0
    return exp(-R/Rd-z/zd);
}

class binney_prior: public galaxy_prior{

    const double thin_norm = 0.0520313; // integral of p(tau) from 0 to 10
    const double thick_norm = 0.15/(13.-8.); // n2/n1 * 1/Delta_tau
    const double halo_factor = pow(R0*R0+z0*z0,3.39/2.);
    const double halo_norm = 0.005 * halo_factor/(13.-10.); // n3/n1 * (R0^2+z0^2)^(3.39/2.) * 1/Delta_tau
    const double thin_metal_mean = 0.;
    const double thin_metal_disp = 0.2;
    const double thin_Rd = 2.6;
    const double thin_zd = 0.3;
    const double thick_metal_mean = -0.6;
    const double thick_metal_disp = 0.5;
    const double thick_Rd = 3.6;
    const double thick_zd = 0.9;
    const double halo_metal_mean = -1.6;
    const double halo_metal_disp = 0.5;
public:
    binney_prior(VecDoub SP):galaxy_prior(SP,true){}
    double prior(VecDoub x, double Z, double tau);
};

class binney_prior_alpha: public galaxy_prior{

    const double thin_norm = 0.0520313; // integral of p(tau) from 0 to 10
    const double thick_norm = 0.15/(13.-8.); // n2/n1 * 1/Delta_tau
    const double halo_factor = pow(R0*R0+z0*z0,3.39/2.);
    const double halo_norm = 0.005 * halo_factor/(13.-10.); // n3/n1 * (R0^2+z0^2)^(3.39/2.) * 1/Delta_tau
    const double thin_metal_mean = 0.;
    const double thin_metal_disp = 0.2;
    const double thin_alpha_mean = 0.;
    const double thin_alpha_disp = 0.1;
    const double thin_Rd = 2.6;
    const double thin_zd = 0.3;
    const double thick_metal_mean = -0.6;
    const double thick_metal_disp = 0.5;
    const double thick_alpha_mean = 0.25;
    const double thick_alpha_disp = 0.1;
    const double thick_Rd = 3.6;
    const double thick_zd = 0.9;
    const double halo_metal_mean = -1.6;
    const double halo_metal_disp = 0.5;
    const double halo_alpha_mean = 0.25;
    const double halo_alpha_disp = 0.25;
public:
    binney_prior_alpha(VecDoub SP):galaxy_prior(SP,true){}
    double prior(VecDoub x, double Z, double tau, double alpha);
    double prior(VecDoub x, double Z, double tau){return prior(x,Z,tau,0.);}
};

class new_prior_2018: public galaxy_prior{

    // From Bland-Hawthorn & Gerhard (2016) (as well as scalelengths below)
    const double local_thick_thin = 0.04;
    const double local_halo_thin = 0.005;

    // Need to normalize age distributions over [0, 12.586 Gyr]
    const double thin_norm = 0.054838; // 1/integral of p(tau)
    const double thick_norm = local_thick_thin*1.108; // n2/n1 * 1/integral of gaussian in age over [0,13]
    const double halo_factor = pow(R0*R0+z0*z0,3.39/2.)*1.624*1.13; // n3/n1 * 1/integral of gaussian in age over [0.,13.]
    // Times 1/integral of p(metal) from -2.2 to 0.6 dex
    const double halo_norm = local_halo_thin * halo_factor;

    const double thin_metal_mean = -0.1;
    const double thin_metal_disp = 0.3;
    const double thin_age_mean = 8.;
    const double thin_age_disp = 1.5;
    const double thin_Rd = 2.6;
    const double thin_zd = 0.3;
    const double thick_metal_mean = -0.6;
    const double thick_metal_disp = 0.5;
    const double thick_tau_mean = 10.;
    const double thick_tau_disp = 2.;
    const double thick_Rd = 2.;
    const double thick_zd = 0.9;
    const double halo_metal_mean = -1.6;
    const double halo_metal_disp = 0.5;
    const double halo_tau_mean = 12.;
    const double halo_tau_disp = 2.;

    const double bulge_age_mean = 10.;
    const double bulge_age_disp = 3.;
    const double bulge_metal_mean = 0.;
    const double bulge_metal_disp = 0.5;

    // Using S model parameters from Simion et al. (2017) with central
    // normalization from Robin et al. (2012), Rc from Sharma et al. (2011)

    const double bulge_phi = 19.57*PI/180.; // in radians
    const double xb0 = 1.47;
    const double yb0 = 0.63;
    const double zb0 = 0.47;
    const double bulge_norm = 35.45/0.04*1.2414*1.135; // central/local thin Msun/pc^3
    // Bovy (2018) total mid-plane stellar density of 0.04 Msun/pc^3
    // Times 1/integral of p(tau) from 0 to 12.586 Gyr
    // Times 1/integral of p(metal) from -2.2 to 0.6 dex
    const double cperp = 1.88;
    const double cpar = 3.04;
    const double bulge_r_core = 2.54;
public:
    new_prior_2018(VecDoub SP):galaxy_prior(SP,true){}
    double thin_prior(double R, double z, double Z, double tau);
    double thick_prior(double R, double z, double Z, double tau);
    double halo_prior(double R, double z, double Z, double tau);
    double bulge_prior(VecDoub x, double Z, double tau);
    double prior(VecDoub x, double Z, double tau);
};
