#include "sfr.h"
//=============================================================================
// Integrands
//=============================================================================
double _gas_consumed_per_unit_radius(double t, void *p){
    gas_con_st *P=(gas_con_st *) p;
    return (*P->SFR)(P->R,t);
}
double _gas_consumed(double R, void *p){
    StarFormationRate *P = (StarFormationRate *) p;
    GaussLegendreIntegrator GL(50.);
    gas_con_st Q = {P,R};
    return 2.*PI*R*GL.integrate(&_gas_consumed_per_unit_radius,0.,P->MaxAge(),&Q);
}
//=============================================================================
// Constructor
StarFormationRate::StarFormationRate(ModelParameters M){
    Rmin=M.parameters["grids"]["MinimumRadius"];
    Rmax=M.parameters["grids"]["MaximumRadius"];
    GalaxyAge=M.parameters["fundamentals"]["GalaxyAge"];
}
//=============================================================================
// Integrations
double StarFormationRate::gas_consumed_per_unit_radius(double R){
    GaussLegendreIntegrator GL(50.);
    gas_con_st Q = {this,R};
    return GL.integrate(&_gas_consumed_per_unit_radius,0.,GalaxyAge,&Q);
}
double StarFormationRate::gas_consumed(void){
    GaussLegendreIntegrator GL(50.);
    return GL.integrate(&_gas_consumed,Rmin,Rmax,this);
}
//=============================================================================
// SFR implementations
double SFR_SB15::operator()(double R, double t){
    double truncation = pow(1.+exp(4.*(R-Rb)),.25*(1./Rd-4./Rd));
    return presentSFR*exp(-t/t_d-t_s/(t+1.e-1))*exp(-R/Rd)*truncation;
}
double SFR_ExpDecay::operator()(double R, double t){
    double truncation = pow(1.+exp(4.*(R-Rb)),.25*(1./Rd-4./Rd));
    return presentSFR*exp(-t/t_d)*exp(-R/Rd)*truncation;
}
double SFR_Neige2020::operator()(double R, double t){
    double truncation = pow(1.+exp(4.*(R-Rb)),.25*(1./Rd-4./Rd));
    double norm = t_sfr/(1-x_io*R/R0)*(exp(-x_io*R/R0*t_m/t_sfr)-exp(-t_m/t_sfr));
    t = GalaxyAge - t + 1e-2;
    double timebit = exp(((1-x_io*R/R0)*t-t_m)/t_sfr)/norm*exp(t_to/GalaxyAge*(t/(t-GalaxyAge)));
    double radiusbit = exp(-R/Rd)*truncation;
    return presentSFR*timebit*radiusbit;
}
//=============================================================================
// Map for creating shared pointer instances of SFR from string of class name
shared_map<StarFormationRate,ModelParameters> sfr_types ={
    {"SB15",&createSharedInstance<StarFormationRate,SFR_SB15>},
    {"Neige2020",&createSharedInstance<StarFormationRate,SFR_Neige2020>},
    {"ExpDecay",&createSharedInstance<StarFormationRate,SFR_ExpDecay>}};
//=============================================================================
