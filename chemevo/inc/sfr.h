#ifndef SFR_H
#define SFR_H
//=============================================================================
#include <map>
#include <string>
//=============================================================================
#include "params.h"
#include "GSLInterface/GSLInterface.h"
#include "utils.h"
//=============================================================================
/**
 * @brief Star Formation rate base class
 */
class StarFormationRate{
private:
protected:
    double Rmin, Rmax; // Minimum and maximum radius of star formation
    double GalaxyAge;
    double presentSFR;
public:
    StarFormationRate(double Rmin, double Rmax, double Tmax)
        :Rmin(Rmin),Rmax(Rmax),GalaxyAge(Tmax),presentSFR(1.){}
    StarFormationRate(ModelParameters M);
    double MaxAge(void){return GalaxyAge;}
    /**
     * @brief Star formation rate
     *
     * @param R radius
     * @param t time
     *
     * @return star formation rate at radius R and time t.
     */
    virtual double operator()(double R, double t)=0;
    /**
     * @brief gas consumed at radius R
     * @details integrated star formation rate at radius R
     *
     * @param R radius
     * @return gas consumed at radius R
     */
    double gas_consumed_per_unit_radius(double R);
    /**
     * @brief gas consumed over all radii and time
     * @details integrated SFR over radius and time
     * @return integrated SFR over radius and time
     */
    double gas_consumed(void);
};

//=============================================================================
/**
 * @brief Star formation rate used in Sanders & Binney (2015)
 */
class SFR_SB15:public StarFormationRate{
private:
    double t_d; // Long time-scale for decline in SFR
    double t_s; // Short time-scale for initial rise in SFR
    double Rd;  // Scale-length of population of stars formed.
    double Rb;  // Radius at which exponential profile is truncated.
public:
    SFR_SB15(double t_d=8.,double t_s=0.5,double Rd=4.,
             double Rmin=0.3, double Rmax=20., double Tmax=12.,double Rb=1000.)
	   :StarFormationRate(Rmin,Rmax,Tmax),t_d(t_d),t_s(t_s),Rd(Rd),Rb(Rb){}
    SFR_SB15(ModelParameters M,double t_d=8.,double t_s=0.5,double Rdx=4.,double Rbx=1000.)
        :StarFormationRate(M),
	 t_d(t_d),
	 t_s(t_s),
	 Rd(Rdx),
     Rb(Rbx){
        Rd=extract_param(M.parameters["fundamentals"],
                         "StarScaleLength", Rdx);
        Rb=extract_param(M.parameters["fundamentals"],
                         "TruncationRadius", Rbx);
        presentSFR=1.;
        presentSFR=(double)M.parameters["fundamentals"]["PresentSFR"]/(
                (*this)(M.parameters["fundamentals"]["SolarRadius"],
                        M.parameters["fundamentals"]["GalaxyAge"]));
        }
    double operator()(double R, double t);
};
class SFR_ExpDecay:public StarFormationRate{
private:
    double t_d; // Long time-scale for decline in SFR
    double Rd;  // Scale-length of population of stars formed.
    double Rb;  // Radius at which exponential profile is truncated.
public:
    SFR_ExpDecay(double t_d=8.,double Rd=4.,
             double Rmin=0.3, double Rmax=20., double Tmax=12.,double Rb=1000.)
       :StarFormationRate(Rmin,Rmax,Tmax),t_d(t_d),Rd(Rd){}
    SFR_ExpDecay(ModelParameters M,double t_dx=8.,double Rdx=4.,double Rbx=1000.)
        :StarFormationRate(M),
     t_d(t_dx),
     Rd(Rdx),
     Rb(Rbx){
        Rd=extract_param(M.parameters["fundamentals"],
                         "StarScaleLength", Rdx);
        Rb=extract_param(M.parameters["fundamentals"],
                         "TruncationRadius", Rbx);
        t_d=extract_param(M.parameters["fundamentals"],
                         "SFR_decay_scale", t_dx);
        presentSFR=1.;
        presentSFR=(double)M.parameters["fundamentals"]["PresentSFR"]/(
                (*this)(M.parameters["fundamentals"]["SolarRadius"],
                        M.parameters["fundamentals"]["GalaxyAge"]));
        }
    double operator()(double R, double t);
};

class SFR_Neige2020: public StarFormationRate{
private:
    double t_sfr;
    double x_io;
    double t_m;
    double t_to;
    double Rd;  // Scale-length of population of stars formed.
    double Rb;  // Radius at which exponential profile is truncated.
    double R0;  // Solar radius.
public:
    SFR_Neige2020(ModelParameters M): StarFormationRate(M), t_m(6.){
        t_sfr = extract_param(M.parameters["fundamentals"],"SFR_decay_scale",1.02);
        x_io = extract_param(M.parameters["fundamentals"],"InsideOutScale",0.67);
        Rd = extract_param(M.parameters["fundamentals"],"StarScaleLength",2.85);
        Rb = extract_param(M.parameters["fundamentals"],"TruncationRadius",1000.);
        t_to = extract_param(M.parameters["fundamentals"],"SFR_ShortTimescale",1.);
	R0 = M.parameters["fundamentals"]["SolarRadius"];
        presentSFR=1.;
        presentSFR=(double)M.parameters["fundamentals"]["PresentSFR"]/(
                (*this)(M.parameters["fundamentals"]["SolarRadius"],
                        M.parameters["fundamentals"]["GalaxyAge"]));
    }
    double operator()(double R, double t);
};

//=============================================================================
// Simple structure for integrating the SFR wrt time.
struct gas_con_st{
    StarFormationRate *SFR;
    double R;
};
//=============================================================================
// Map for creating new instances of SFR from a string giving the class name
// sfr_types[std::string s="<ClassName>"](ModelParameters M)
// produces a shared pointer of Class <ClassName> initialized with the model
// parameters M.
extern shared_map<StarFormationRate,ModelParameters> sfr_types;
//=============================================================================
#endif
//=============================================================================
