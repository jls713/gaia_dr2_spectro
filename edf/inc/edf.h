#ifndef EDF_
#define EDF_

// ===========================================================================
#include "utils.h"
#include "json.hpp"
using json = nlohmann::json;
#include "cuba/cuba.h"
#include "potential.h"
#include "tables_aa.h"
#include "edf_base.h"
#include "metal_relation.h"
#include "rad_mig.h"
#include "sfr.h"
#include "thin_disc.h"
#include "thick_disc.h"
#include "halo.h"
// ===========================================================================
// In Sanders & Binney (2015), the EDF is chosen to describe the distribution
// function of all stars that have ever existed -- therefore, Gamma(tau) is
// truly the star formation history. When multiplying by the IMF, we produce
// a full model for all stars that have ever existed. Our selection function
// then, in addition to computing the relative fraction of different stellar
// populations, also describes the observable fraction today (i.e. those still
// alive).
//
// There is a slight inconsistency where the EDF without multiplying by the
// selection function now includes all the stars that have ever lived, and so
// incorrectly includes a fraction of stars that would have died. This
// fraction obviously increases with age. However, the fraction of a stellar
// population that has died is quite small and we assume the stellar
// distribution is dominated by the low mass stars. This was, however, an
// oversight.
//
// A different formalism would be for the EDF to describe the current stellar
// distribution. Then it is correct to compute the density etc. today from
// this. The selection fraction would then be the integral of the IMF over the
// selection divided by the integral of the IMF over surviving stars (or
// rather, the integral over the present mass function, PMF). In this
// case, Gamma(tau) is no longer the star formation history, but the present
// day age distribution.
//
// For further modelling, it is easier to work with the former definition. The
// EDF is then multipled by the IMF without any normalization that is age (or
// metallicity dependent). In this way, we derive the star formation history
// which can be used in further chemical evolution modelling.
// ===========================================================================

class sb15_edf: public edf{
private:
	RadialMigration *rm;
	metal_relation *mr;
	double totalDiscMass;
	SFR_base *sfr;
	VecDoub StandardSolar;

	double tau_m=12.,tau_1=0.11,T0=8.,tau_T=10.;
    double BETA_R=0.33, BETA_Z=0.4;
	double SIGHALO=1.,FHALO=-1.5,SIGFHALO=0.5, GAMMAT=0.5;

public:
	thin_disc_edf *thind;
	thick_disc_edf *thickd;
	halo_edf *halod;
	bool ifthin, ifthick, ifhalo;

	sb15_edf(Potential_JS *Pot, Actions_AxisymmetricFudge_InterpTables *Act, VecDoub params= {3.6,9.33,0.36,0.25,3.,5.3,0.4,0.58,0.2,5.7,2.4,-1.,7.5,-0.08,1.,0.},VecDoub StandardSolar=conv::StandardSolarPAUL, bool NewMetalRelation=true,bool ContinuousSFR=true,int order=20);
	inline std::string name(void){return "Sanders-Binney (2015) EDF";}
	inline VecDoub SunCoords(void){return StandardSolar;}


	// Evaluate full EDF functions
	double fullDF_actions(const VecDoub& Actions, double age, double RcP);
	double fullDF_actions_Z(const VecDoub& Actions, double age, double Z);
	double chemDF_actions(const VecDoub& Actions, double F);
	double chemDF_real(const VecDoub& X, double F);
	double chemDF_actions_justdiscs(const VecDoub& Actions, double FeH);
	double chemDF_real_justdiscs(const VecDoub& X, double FeH);
	double full_DF(const VecDoub& X, double age, double RcP);
	double full_DF_Z(const VecDoub& X, double age, double Z);

	double normalization(double scale=1.);
	inline void reset_pot(Potential_JS *Pot){pot=Pot;ActionCalculator->partial_reset(pot);}
	inline void setdiscmass(double dm){totalDiscMass=dm;sfr->reset_mass(dm);}

	// Manipulate parameters
	std::string params(void);
	void readParams(std::string);
	void setParams(VecDoub A, bool print=false);
	void setParams(std::string jsonfile, bool print=false);
	void printParams(std::ostream&out=std::cout);
	void simpleprintParams(std::ofstream *output_file=nullptr);
	VecDoub returnParams(void);

	// Turn On and Off components
	inline void TurnOffThin(){ifthin = 0;}
	inline void TurnOnThin(){ifthin = 1;}
	inline void TurnOffThick(){ifthick = 0;}
	inline void TurnOnThick(){ifthick = 1;}
	inline void TurnOffHalo(){ifhalo = 0;}
	inline void TurnOnHalo(){ifhalo = 1;}

	// Metallicity and age properties
	inline double get_taum(void){return sfr->tau_m;}
	inline double FeH(double a, double b){return mr->FeH(a,b);}
	inline double maxZ(void){return mr->FeH(0.,0.);}
	inline double minZ(void){return mr->FEHStart;}
	inline double AgeFromMetalRadius(double a, double b){return mr->AgeFromMetalRadius(a,b);}
	inline double RadiusFromMetal(double a, double b){return mr->RadiusFromMetal(a,b);}
	inline double MaxAge(double Z){return mr->MaxAge(Z);}
	inline double get_min_Z(void){
		if(ifhalo)
			return halod->FHALO-4.*halod->SIGFHALO;
		else
			return minZ();
	}
	inline double get_max_Z(void){
		if(ifhalo and halod->FHALO+4.*halod->SIGFHALO>maxZ())
			return halod->FHALO+4.*halod->SIGFHALO;
		else
			return maxZ();
	}

	// Moments
	double density(double R, double z, double IE, VecDoub FF = {-10.,10.});
	double FHist(double R, double z, double F, double IE);
	double UHist(double R, double z, double v, double IE);
	double VHist(double R, double z, double v, double IE);
	double WHist(double R, double z, double v, double IE);
	VecDoub mean_proper_motions(double l, double b, double s, double IE=1e-3);
	double integrate_along_line_of_sight(double R, double z, double IE);
	VecDoub moments(double R, double z, double tau, double Z, double IE=1e-3);
};

struct edf_norm_struct{
	sb15_edf *edf;
	VecDoub x2min,x2max;
	double scale;
	edf_norm_struct(sb15_edf *edf, VecDoub x2min, VecDoub x2max,double scale)
		:edf(edf),x2min(x2min),x2max(x2max),scale(scale){}
};

struct edf_density_struct: edf_norm_struct{
	double R, Z;
	edf_density_struct(sb15_edf *edf, double r, double z, VecDoub x2m, VecDoub x2n)
		:edf_norm_struct(edf,x2m,x2n,1.),R(r),Z(z){}
};

struct edf_vhist_struct: edf_density_struct{
	double V;
	edf_vhist_struct(sb15_edf *edf, double r, double z, double v, VecDoub x2m, VecDoub x2n)
		:edf_density_struct(edf,r,z,x2m,x2n),V(v){}
};

struct edf_moments_struct: edf_density_struct{
	double met, tau;
	edf_moments_struct(sb15_edf *edf, double R, double z, double met, double tau, VecDoub x2m, VecDoub x2n)
		:edf_density_struct(edf,R,z,x2m,x2n),met(met),tau(tau){}
};

double edf_integrate(integrand_t integrand, edf_norm_struct *P, double IE, double AE=0, std::string type="Divonne", double *err=nullptr, int verbosity=0, std::vector<VecDoub> *peaks=nullptr);
VecDoub edf_multi_integrate(integrand_t integrand, int nintegrals, edf_norm_struct *P, double IE, double AE=0, std::string type="Divonne", double *err=nullptr, int verbosity=0);

json edf_parameters(void);

#endif

// ===========================================================================
