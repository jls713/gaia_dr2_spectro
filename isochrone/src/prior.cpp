#include "prior.h"
//=============================================================================

double binney_prior::prior(VecDoub X, double Z, double tau){
    double R = sqrt(X[0]*X[0]+X[1]*X[1]), z = X[2];
    if(tau>13.) return 0.;
    double T = 0.;
    if(tau<10.) T+= thin_norm*
                    GFunction(Z-thin_metal_mean,thin_metal_disp)*
                    exp(0.119*tau)*
                    double_exponential(R-R0,fabs(z)-z0,thin_Rd,thin_zd);
    if(tau>8.) T+= thick_norm*
                    GFunction(Z-thick_metal_mean,thick_metal_disp)*
                    double_exponential(R-R0,fabs(z)-z0,thick_Rd,thick_zd);
    if(tau>10.) T+= halo_norm*
                    GFunction(Z-halo_metal_mean,halo_metal_disp)*
                    pow(R*R+z*z,-3.39/2.);
    return T;
}
//=============================================================================

double binney_prior_alpha::prior(VecDoub X, double Z, double tau, double alpha){
    double R = sqrt(X[0]*X[0]+X[1]*X[1]), z = X[2];
    if(tau>13.) return 0.;
    double T = 0.;
    if(tau<10.) T+= thin_norm*
                    GFunction(Z-thin_metal_mean,thin_metal_disp)*
                    GFunction(alpha-thin_alpha_mean,thin_alpha_disp)*
                    exp(0.119*tau)*
                    double_exponential(R-R0,fabs(z)-z0,thin_Rd,thin_zd);
    if(tau>8.) T+= thick_norm*
                    GFunction(Z-thick_metal_mean,thick_metal_disp)*
                    GFunction(alpha-thick_alpha_mean,thick_alpha_disp)*
                    double_exponential(R-R0,fabs(z)-z0,thick_Rd,thick_zd);
    if(tau>10.) T+= halo_norm*
                    GFunction(Z-halo_metal_mean,halo_metal_disp)*
                    GFunction(alpha-halo_alpha_mean,halo_alpha_disp)*
                    pow(R*R+z*z,-3.39/2.);
    return T;
}
//=============================================================================

double sech2(double x){
    return pow(2./(exp(x)+exp(-x)),2.);
}
double core_function(double x, double y, double rc){
    double r = sqrt(x*x+y*y);
    if(r<rc) return 1.;
    else return exp(-pow((r-rc)/0.5,2.));
}

double new_prior_2018::bulge_prior(VecDoub X, double Z, double tau){
    double xrot = X[0]*cos(bulge_phi)-X[1]*sin(bulge_phi);
    double yrot = -X[0]*sin(bulge_phi)-X[1]*cos(bulge_phi);
    return bulge_norm*
           GFunction(Z-bulge_metal_mean,bulge_metal_disp)*
           GFunction(tau-bulge_age_mean,bulge_age_disp)*
 	       sech2(pow(pow(fabs(xrot/xb0),cperp)+
                     pow(fabs(yrot/yb0),cperp),cpar/cperp)+
                 pow(fabs(X[2]/zb0),cpar))*
           core_function(X[0],X[1],bulge_r_core);
}
double new_prior_2018::thin_prior(double R, double z, double Z, double tau){
    double thin_age = 0.;
    if(tau<thin_age_mean)
        thin_age = exp(0.119*tau);
    else
        thin_age = exp(-pow((tau-thin_age_mean)/thin_age_disp,2.)*.5)*exp(0.119*thin_age_mean);
    return thin_norm*
           GFunction(Z-thin_metal_mean,thin_metal_disp)*
           thin_age*
           double_exponential(R-R0,fabs(z)-z0,thin_Rd,thin_zd);
}
double new_prior_2018::thick_prior(double R, double z, double Z, double tau){
    return thick_norm*
           GFunction(Z-thick_metal_mean,thick_metal_disp)*
           GFunction(tau-thick_tau_mean,thick_tau_disp)*
           double_exponential(R-R0,fabs(z)-z0,thick_Rd,thick_zd);
}
double new_prior_2018::halo_prior(double R, double z, double Z, double tau){
    return halo_norm*
           GFunction(Z-halo_metal_mean,halo_metal_disp)*
           GFunction(tau-halo_tau_mean,halo_tau_disp)*
           pow(R*R+z*z+1e-2,-3.39/2.);
}
double new_prior_2018::prior(VecDoub X, double Z, double tau){
    double R = sqrt(X[0]*X[0]+X[1]*X[1]), z = X[2];
    if(tau>13.) return 0.;
    return thin_prior(R,z,Z,tau)+thick_prior(R,z,Z,tau)+halo_prior(R,z,Z,tau)+bulge_prior(X,Z,tau);
}
//=============================================================================
