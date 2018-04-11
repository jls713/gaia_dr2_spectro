#include "imf.h"
//=============================================================================
double _imf(double m, void *p){
    InitialMassFunction *P = (InitialMassFunction *) p;
    return (*P)(m);
}

double InitialMassFunction::norm(void){
    GaussLegendreIntegrator GL(50.);
    return GL.integrate(&_imf,MinimumMass,MaximumMass,this);
}

double _imf_m(double m, void *p){
    InitialMassFunction *P = (InitialMassFunction *) p;
    return m*(*P)(m);
}

InitialMassFunction::InitialMassFunction(ModelParameters M){
    auto F = M.parameters["fundamentals"];
    std::vector<std::string> vars = {"MinimumMass","MaximumMass"};
    for(auto v:vars)
        if (F.find(v) == F.end()) {
            LOG(INFO)<<v<<" not found in parameters file\n";
            throw std::invalid_argument(v+" not found in parameters file");
        }
    MinimumMass = M.parameters["fundamentals"]["MinimumMass"];
    MaximumMass = M.parameters["fundamentals"]["MaximumMass"];
}

double InitialMassFunction::mean_mass(void){
    GaussLegendreIntegrator GL(50.);
    return GL.integrate(&_imf_m,MinimumMass,MaximumMass,this);
}
//=============================================================================
double SalpeterIMF::operator()(double m){
    return pow(m,-2.35)/Norm;
}
//=============================================================================
double TinsleyIMF::operator()(double m){
    if(m<2.) return pow(m,-2.)/Norm;
    else if(m<10.) return pow(m,-2.3)/Norm;
    else return 10.*pow(m,-3.3)/Norm;
}
//=============================================================================
double ScaloIMF::operator()(double m){
    if(m<1.) return 0.39*pow(m,-1.2)/Norm;
    else if(m<10.) return 0.39*pow(m,-2.7)/Norm;
    else return 0.1553*pow(m,-2.3)/Norm;
}
//=============================================================================
double KroupaIMF::operator()(double m){
    if(m<0.08) return 0.4375*pow(m,-0.3)/Norm;
    else if(m<0.5) return 0.035*pow(m,-1.3)/Norm;
    else if(m<1.) return 0.019*pow(m,-2.3)/Norm;
    else return 0.019*pow(m,-2.7)/Norm;
}
//=============================================================================
double ChabrierIMF::operator()(double m){
    if(m<1.) return exp(-pow(log(m/0.2),2.)/0.6)/Norm;
    else return pow(m,-2.35)/Norm;
}
//=============================================================================
// Map for creating new instances of IMF from string.
shared_map<InitialMassFunction,ModelParameters> imf_types ={
    {"Salpeter",&createSharedInstance<InitialMassFunction,SalpeterIMF>},
    {"Tinsley",&createSharedInstance<InitialMassFunction,TinsleyIMF>},
    {"Scalo",&createSharedInstance<InitialMassFunction,ScaloIMF>},
    {"Kroupa",&createSharedInstance<InitialMassFunction,KroupaIMF>},
    {"Chabrier",&createSharedInstance<InitialMassFunction,ChabrierIMF>}
};
//=============================================================================
