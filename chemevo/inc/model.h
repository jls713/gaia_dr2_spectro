#ifndef MODEL_H
#define MODEL_H
//=============================================================================
#include <utility>
#include <map>
#include <memory>
#include <exception>
//=============================================================================
#include <omp.h>
#include "easylogging++.h"
#include <H5Cpp.h>
//=============================================================================
#include "utils_ch.h"
#include "params.h"
#include "sfr.h"
#include "imf.h"
#include "ages.h"
#include "iarates.h"
#include "yields.h"
#include "in_out.h"
#include "grid.h"
#include "radmig.h"
//=============================================================================
class Model{
private:
	//=========================================================================
	ModelParameters params;
	std::shared_ptr<SolarAbundances> solar;
	//=========================================================================
	std::shared_ptr<StarFormationRate> sfr;
	std::shared_ptr<InitialMassFunction> imf;
	std::shared_ptr<StellarLifetime> lifetime;
	std::shared_ptr<RadialMigration> rad_mig;
	//=========================================================================
	std::unique_ptr<TypeIaRate> typeIaRate;
	std::unique_ptr<YieldsSet> yields;
	std::unique_ptr<Inflow> inflow;
	std::unique_ptr<Outflow> outflow;
	std::unique_ptr<RadialFlow> radialflow;
	std::unique_ptr<GasDump> gasdumpflow;
	//=========================================================================
	std::unique_ptr<Grid> gas_mass; // total gas mass
	std::unique_ptr<Grid> warm_gas_mass; // total warm gas mass
	std::unique_ptr<Grid> stellar_mass; // total stellar mass
	std::unique_ptr<Grid> reducedSFR; // SFR-DeathRate
	std::vector<Grid> mass_fraction; // First two elements H and He, then all other required elements
	std::vector<Grid> mass_fraction_warm; // First two elements H and He, then all other required elements
	std::unique_ptr<Grid> metallicity; // Z = 1-X-Y = 1-first two entries of mass_fraction
	//=========================================================================
	std::map<int,Element> elements;
	std::map<Element,int> elements_r;
	//=========================================================================
	bool agb_yields, typeII_yields, typeIa_yields, migration;
    bool single_zone, gasdump=false, use_warm_phase, logspace;
    double warm_cold_ratio, warm_cooling_time;
    //=========================================================================
	unsigned iteratemax; // maximum number of iterations to perform per step
	double tol;          // target relative accuracy between iterations
protected:
public:
	// Constructors
	Model(std::string params_file):params(params_file){
		setup();
		fill_initial_grids();
	}
	Model(ModelParameters &params):params(params){
		setup();
		fill_initial_grids();
	}
	//=========================================================================
	// Functions for setting up and running models
	void setup(void);
	void fill_initial_grids(void);
	int check_parameters(void);
	int step(unsigned nt, double dt);
	void simple_step(unsigned nt, double dt);
	void run(void);
	void write(std::string filename);
	void write_properties(H5::H5File &fout);
	void expand_grids(unsigned nt, double t);
        int check_metallicity(double R, double t, double dt, unsigned nR, unsigned NR, unsigned nt);
	//=========================================================================
	ModelParameters &parameters(void){return params;}
	//=========================================================================
	double X(double R, double t);
	double Y(double R, double t);
	double Z(double R, double t);
	//=========================================================================
	//=========================================================================
	// Supernovae rates
	double SNIaRate(double R, double t);
	double SNIIRate(double R, double t);
	//=========================================================================
	// Number of stars dying per unit time
	double DeathRate(double R, double t);
	double SwitchDeathRate(double R, double t, double tmin, double tmax, std::string s);
	double TestSNIIRate(double R, double t);
	double AGBDeathRate(double R, double t);
	//=========================================================================
	// Enrichment rates of element E at radius R and time t
	double EnrichmentRate(Element E, double R, double t);
	double TypeIIEnrichmentRate(Element E, double R, double t);
	double TypeIaEnrichmentRate(Element E, double R, double t);
	double AGBEnrichmentRate(Element E, double R, double t);
	double SwitchEnrichmentRate(Element E, double R, double t, double tmin, double tmax, std::string s);
	//=========================================================================
	// Mass of gas returned at radius R and time t
	double GasReturnRate(double R, double t);
	double TypeIIGasReturnRate(double R, double t);
	double TypeIaGasReturnRate(double R, double t);
	double AGBGasReturnRate(double R, double t);
	double SwitchGasReturnRate(double R, double t, double tmin, double tmax, std::string s);
	//=========================================================================
	// Interfaces to member classes
	double SFR(double R, double t);
	double IMF(double m);
	double Lifetime(double m,double Z);
	double EjectedMass(Element E, double M, double Z);
	double TotalEjectedMass(double M, double Z);
	double Mass_Star_Dying_Now(double age,double Z);
	double dMdt(double M, double age, double Z);
	double InflowRate(double R, double t);
	double OutflowRate(double R, double t, double SFR, double GasReturn);
	double RadialFlowRate(double m, double mup, double R, double Rdown,
	                      double Rup, double t, double dt, int*err);
	double RadialFlowRateFromGrid(double R, double t, double dt,
	                              double gm_prev, unsigned nR, unsigned nt,
	                              unsigned NR, int*err);
	double EnrichRadialFlowRateFromGrid(Element E, double R, double t,
	                                    double dt, double e_gm_prev,
	                                    unsigned nR, unsigned nt,
	                                    unsigned NR, int*err);
	double RadialMigrationKernel(double R, double Rp, double t);
	double GasDumpRate(double R, double t, double dt);
	double EnrichGasDumpRate(Element E, double R, double t, double dt);
	//=========================================================================
	// Printing functions
	void print_abundance_grid(std::ostream &out,Element E);
	void print_abundance_grid_fixed_radius(std::ostream &out,Element E,double R);
	void print_abundance_grid_fixed_time(std::ostream &out,Element E,double t);
	void print_abundance_grid_now(std::ostream &out,Element E);

	void print_gasmass_grid(std::ostream &out);
	void print_gasmass_grid_fixed_radius(std::ostream &out,double R);
	void print_gasmass_grid_fixed_time(std::ostream &out,double t);
	void print_gasmass_grid_now(std::ostream &out);

	void print_metallicity_grid(std::ostream &out);
	void print_metallicity_grid_fixed_radius(std::ostream &out,double R);
	void print_metallicity_grid_fixed_time(std::ostream &out,double t);
	void print_metallicity_grid_now(std::ostream &out);
	//=========================================================================
};
//=============================================================================
// Simple structure for computing integrals of rates
struct Rate_st{
    Model *M;
    double R;
    double t;
    std::string which;
    double mthresh;
    Rate_st(Model *M, double R, double t, std::string which, double mthresh):
    	M(M),R(R),t(t),which(which),mthresh(mthresh){}
};
// Simple structure for computing integrals of rates with radial migration
struct Rate_st_2D:Rate_st{
    VecDoub x2min,x2max;
    Rate_st_2D(Model *M, double R, double t, VecDoub x2min, VecDoub x2max,
    std::string which, double mthresh)
	:Rate_st(M,R,t,which,mthresh),x2min(x2min),x2max(x2max){}
};
// Simple structure for computing enrichment rate integrals
struct EnrichRate_st:Rate_st{
    Element E;
    EnrichRate_st(Model *M,double R,double t,Element E,
    std::string which, double mthresh):Rate_st(M,R,t,which,mthresh),E(E){}
};
// Simple structure for computing enrichment rate integrals with radial mig.
struct EnrichRate_st_2D:EnrichRate_st{
    VecDoub x2min,x2max;
    EnrichRate_st_2D(Model *M, double R, double t, Element E, VecDoub x2min, VecDoub x2max, std::string which,double mthresh)
	:EnrichRate_st(M,R,t,E,which,mthresh),x2min(x2min),x2max(x2max){}
};
//=============================================================================
#endif
//=============================================================================
