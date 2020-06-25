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
double truncfunc(double R, double Rd, double Rtrunc){
    return pow(1.+exp(4.*(R-Rtrunc)),.25*(1./Rd-4./Rd));
}
double SFR_SB15::operator()(double R, double t){
    return presentSFR*exp(-t/t_d-t_s/(t+1.e-1))*exp(-R/Rd)*truncfunc(R, Rd, Rb);
}
double SFR_ExpDecay::operator()(double R, double t){
    return presentSFR*exp(-t/t_d)*exp(-R/Rd)*truncfunc(R, Rd, Rb);
}
double SFR_DoubleInfall::operator()(double R, double t){
    double sfr = original_SFR_scaling*exp(-t/t_d)*exp(-R/Rd)*truncfunc(R, Rd, Rb);
    //if(t>=gasdumpflow->dump_time()) sfr+=coeff*KSA*pow((*gasdumpflow)(R,gasdumpflow->dump_time(),1.),KSN)*exp(-(t-gasdumpflow->dump_time())/t_d_2);
    sfr+=coeff*KSA*pow((*gasdumpflow)(R,t,1.),KSN)*exp(-(t-gasdumpflow->dump_time())/t_d_2);
    return presentSFR*sfr;
}
double SFR_Neige2020::operator()(double R, double t){
    double norm = t_sfr/(1-x_io*R/R0)*(exp(-x_io*R/R0*t_m/t_sfr)-exp(-t_m/t_sfr));
    t = GalaxyAge - t + 1e-2;
    double timebit = exp(((1-x_io*R/R0)*t-t_m)/t_sfr)/norm*exp(t_to/GalaxyAge*(t/(t-GalaxyAge)));
    double radiusbit = exp(-R/Rd)*truncfunc(R, Rd, Rb);
    return presentSFR*timebit*radiusbit;
}
//=============================================================================
// Map for creating shared pointer instances of SFR from string of class name
shared_map<StarFormationRate,ModelParameters> sfr_types ={
    {"SB15",&createSharedInstance<StarFormationRate,SFR_SB15>},
    {"DoubleInfall",&createSharedInstance<StarFormationRate,SFR_DoubleInfall>},
    {"Neige2020",&createSharedInstance<StarFormationRate,SFR_Neige2020>},
    {"ExpDecay",&createSharedInstance<StarFormationRate,SFR_ExpDecay>}};
//=============================================================================
