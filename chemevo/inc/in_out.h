#ifndef IN_OUT_H
#define IN_OUT_H
//=============================================================================
#include "params.h"
#include "solar.h"
//#include "sfr.h"
#include "grid.h"
//=============================================================================
/**
 * @brief Base class for implementing inflows of gas
 */
class Inflow{
protected:
	VecDoub mass_fraction; // store abundances of inflowing gas
	std::map<Element, int> elements;
	std::shared_ptr<SolarAbundances> solar;
public:
	/**
	 * @brief inflow rate at radius R and time t
	 *
	 * @param R radius
	 * @param t time t
	 *
	 * @return inflow rate at radius R and time t
	 */
	Inflow(ModelParameters M, std::shared_ptr<SolarAbundances> solar,
	       double present_rSFR=0.);
	virtual double operator()(double R, double t, Grid *rSFR=nullptr)=0;
	double Xi_inflow(Element e) {return mass_fraction[elements[e]];};
};
/**
 * @brief No inflow model
**/
class InflowNone: public Inflow{
public:
	InflowNone(ModelParameters M,
	           std::shared_ptr<SolarAbundances> solar,
		       double present_rSFR=0.)
	: Inflow(M, solar, present_rSFR){}
	double operator()(double R, double t, Grid *rSFR=nullptr){return 0.;}
};
/**
 * @brief Simple double inflow model
 * @details gas infall at two rates -- one fast initial rate and a slower more
 * continuous infall. The gas falls in according to an exponential profile
 *
 */
class DoubleInfallInflow: public Inflow{
	double tf; // fast time-scale
	double ts; // slow time-scale
	double Rd; // gas scale-length
	double weight;// relative weight
	double pif; // present infall rate
public:
	DoubleInfallInflow(ModelParameters M,
	                   std::shared_ptr<SolarAbundances> solar,
	                   double present_rSFR=0.);
	double operator()(double R, double t, Grid *rSFR=nullptr);
};
//=============================================================================
/**
 * @brief Base class for implementing outflows of produced gas
 */
class Outflow{
public:
	/**
	 * @brief outflow at radius R and time t
	 *
	 * @param R radius
	 * @param t time t
	 *
	 * @return outflow fraction at radius R and time t
	 */
	virtual double operator()(double R, double t, double SFR, double GasReturn)=0;
};
/**
 * @brief No outflow model
**/
class OutflowNone: public Outflow{
public:
	OutflowNone(ModelParameters M, double prSFR=0.){}
	double operator()(double R, double t, double SFR, double GasReturn){return 0.;}
};
/**
 * @brief Simple Galactic Fountain
 * @details simple outflow model that has two constant outflow fractions (one
 * inside the radius transR and one outside).
 */
class SimpleGalacticFountain: public Outflow{
private:
	double feject_in; // Fraction of products ejected for R<transR
	double feject_out;// Fraction of products ejected for R>transR
	double transR;    // Transition radius
	double dR;        // Radial range over which change occurs
	double A, B; 	  // = (1/2)(fout+fin), (1/2)(fout-fin)
public:
	SimpleGalacticFountain(ModelParameters M);
	double operator()(double R, double t, double SFR, double GasReturn);
};
/**
 * @brief Entrained wind
 * @details simple outflow model for cold gas: outflow = eta_wind SFR * exp(-(R-Rsc)/Rsc)
 */
class EntrainedWind: public Outflow{
private:
	double eta_wind, scale_length;
public:
	EntrainedWind(ModelParameters M);
	double operator()(double R, double t, double SFR, double GasReturn);
};

//=============================================================================
/**
 * @brief Base class for implementing radial flows
 */
class RadialFlow{
public:
	RadialFlow(ModelParameters M,
               std::shared_ptr<SolarAbundances> solar=nullptr,
               double pSFR=0.){}
	/**
	 * @brief flow rate of gas at radius R and time t
	 *
	 * @param R radius
	 * @param t time
	 *
	 * @return flow rate of gas at radius R and time t
	 */
	virtual double flow_rate(double R, double t, Grid *rSFR=nullptr)=0;
	/**
	 * @brief beta_g -- mass loss rate per unit mass in annulus R
	 *
	 * @param R radius R
	 * @param Rdown radius of inner annulus
	 * @param Rup  radius of outer annulus
	 * @param t time t
	 * @param dt timestep
	 * @param err returns 1 if Courant condition not satisfied
	 * @return mass loss rate per unit mass in annulus R
	 */
	double beta_g(double R, double Rdown, double Rup, double t, double dt, Grid *rSFR=nullptr, int*err=nullptr);
	/**
	 * @brief gamma_g -- mass gain rate per unit mass in annulus R from upper
	 * annulus
	 *
	 * @param R radius R
	 * @param Rdown radius of inner annulus
	 * @param Rup  radius of outer annulus
	 * @param t time t
	 * @param dt timestep
	 * @param err returns 1 if Courant condition not satisfied
	 * @return mass loss rate per unit mass in annulus R
	 */
	double gamma_g(double R, double Rdown, double Rup, double t, double dt, Grid *rSFR=nullptr, int*err=nullptr);
	/**
	 * @brief rate of change of mass in annulus R
	 *
	 * @param mass mass in annulus R
	 * @param massup mass in outer annulus
	 * @param R radius
	 * @param Rdown radius of inner annulus
	 * @param Rup radius of outer annulus
	 * @param t time t
	 * @param dt timestep
	 * @param err returns 1 if Courant condition not satisfied
	 * @return rate of change of mass in annulus R
	 */
	double dMdt(double mass, double massup, double R, double Rdown, double Rup, double t, double dt, Grid *rSFR=nullptr, int*err=nullptr);
};
/**
 * @brief Simple linear radial flow
 */
class RadialFlowNone:public RadialFlow{
private:
public:
	RadialFlowNone(ModelParameters M,
	               std::shared_ptr<SolarAbundances> solar=nullptr,
	               double present_rSFR=0.)
		:RadialFlow(M,solar,present_rSFR){}
	double flow_rate(double R, double t, Grid *rSFR=nullptr){
		return 0.;
	}
	double dMdt(double mass, double massup, double R, double Rdown, double Rup, double t, double dt, Grid *rSFR=nullptr, int*err=nullptr){
		return 0.;
	}
	double beta_g(double R, double Rdown, double Rup, double t, double dt, Grid *rSFR=nullptr, int*err=nullptr){return 0.;}
	double gamma_g(double R, double Rdown, double Rup, double t, double dt, Grid *rSFR=nullptr, int*err=nullptr){return 0.;}
};
/**
 * @brief Simple linear radial flow
 */
class LinearRadialFlow:public RadialFlow{
private:
	double VGrad; // Gradient of gas flow velocity wrt radius
public:
	LinearRadialFlow(ModelParameters M,
	                 std::shared_ptr<SolarAbundances> solar=nullptr,
	                 double present_rSFR=0.);
	double flow_rate(double R, double t, Grid *rSFR=nullptr){
		return -VGrad*R;
	}
};
/**
 * @brief Pezzulli radial flow reduced SFR -- templated such that it can be
 * used as either an inflow or a radial flow (as both use the same code)
 */
template<class T>
class PezzulliInflowRadialFlow_rSFR:public T{
private:
	double Alpha; // Dimensionless angular momentum lag of infalling material w.r.t. disc
	double AlphaGrad; // Gradient of the angular momentum lag with radius
	double KSN;   // Kennicutt-Schmidt power
	double A;     // Kennicutt-Schmidt coefficient -- sfr = A \sigma_g^KSN
public:
	PezzulliInflowRadialFlow_rSFR(ModelParameters M,
	                              std::shared_ptr<SolarAbundances> solar,
	                              double present_rSFR=0.);
	~PezzulliInflowRadialFlow_rSFR(){}
	double KSCoeff(double sfr);
	double flow_rate(double R, double t, Grid *rSFR);
	double radial_flow_rate(double R, double t, Grid *rSFR);
	double alpha(){return Alpha;}
	double alphagrad(){return AlphaGrad;}
	double alpha_fn(double R);
	double sigmagas(double R, double t, Grid *rSFR);
	double sigmaeffdot(double R, double t, Grid *rSFR);
	double acc_rate(double R, double t, Grid *rSFR);
	double mu(double R, double t, Grid *rSFR);
	double dmudR(double R, double t, Grid *rSFR);
	double operator()(double R, double t, Grid *rSFR);
};
template<class T>
struct mu_rSFR_st{
    PezzulliInflowRadialFlow_rSFR<T> *P;
    double t;
    Grid *rSFR;
};
//=============================================================================
class GasDump: public Inflow{
public:
	/**
	 * @brief gas dump at radius R and time t
	 *
	 * @param R radius
	 * @param t time t
	 *
	 * @return surface density of gas dumped at radius R, time t
	 */
	GasDump(ModelParameters M, std::shared_ptr<SolarAbundances> solar)
	:Inflow(M, solar){}
	// fudge so we can use inflow class
	double operator()(double R, double t, Grid* rSFR=nullptr){return 0.;}
	virtual double operator()(double R, double t, double dt)=0;
	virtual double dump_time(void)=0;
	// virtual double elements(Element E, double R, double t, double dt) = 0;
};
class GasDumpNone: public GasDump{
public:
	GasDumpNone(ModelParameters M, std::shared_ptr<SolarAbundances> solar)
	: GasDump(M,solar){}
	double operator()(double R, double t, double dt){ return 0.;}
	double dump_time(void){return 0.;}
	// double elements(Element E, double R, double t, double dt){ return 0.;}
};
class GasDumpSimple: public GasDump{
private:
	double surfacedensity, time, central_radius, radial_width;
	double metallicity, alpha;
public:
	GasDumpSimple(ModelParameters M, std::shared_ptr<SolarAbundances> solar);
	double dump_time(void){return time;}
	double operator()(double R, double t, double dt);
	// double elements(Element E, double R, double t, double dt);
};
//=============================================================================
// Maps for creating unique pointers to flow rate classes using string of
// class name
extern unique_map<Inflow,ModelParameters,
            std::shared_ptr<SolarAbundances>,double> inflow_types;
extern unique_map<Outflow,ModelParameters> outflow_types;
extern unique_map<RadialFlow,ModelParameters,
            std::shared_ptr<SolarAbundances>,double> radialflow_types;
extern unique_map<GasDump,ModelParameters,
            std::shared_ptr<SolarAbundances>> gasdump_types;
//=============================================================================
#endif


