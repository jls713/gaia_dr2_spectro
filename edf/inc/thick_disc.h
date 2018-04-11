#ifndef THICK_DISCH
#define THICK_DISCH
// ===========================================================================
// Thick EDF
// ===========================================================================
#include "edf_base.h"
#include "sfr.h"
#include "metal_relation.h"
// ===========================================================================
class thick_disc_edf: public edf{
protected:
	metal_relation *mr;
	SFR_base *sfr;
	RadialMigration *rm;
	const bool usecombinationLzJzforSigma = false;
	bool use_average_radius = false;
	const double L0 = 10.;
	GaussLegendreIntegrator GLI;
	VecDoub StandardSolar;
	double Norm;
public:
	double RD,RSigma,sigmaR,sigmaZ,tau_T,tau_m,RSigmaZ;
	double RD_IO=0., sigmaR0=0., sigmaZ0=0., weight=0.;
	bool use_original_radius_for_drift = false; // if true, use the radius of the population in the drift term for the radial migration -- preserves the disc profile
	bool local_norm=false; // if true, the SFR is the surface density SFR at the Sun (averaged over the column), else it is the global star formation rate per unit mass
	bool continuous_sfr=true; // if true, use continuous SFR Gamma(tau) else use the weight
	thick_disc_edf(Potential_JS *Pot, Actions_AxisymmetricFudge_InterpTables *Act,metal_relation *mr,SFR_base *sfr,RadialMigration *rm, VecDoub params,int order=10):edf(Pot,Act),mr(mr),sfr(sfr),rm(rm),GLI(order),RD(params[0]),RSigma(params[1]),sigmaR(params[2]),sigmaZ(params[3]),tau_T(params[4]),tau_m(params[5]),RSigmaZ(0.),RD_IO(0.),sigmaR0(0.),sigmaZ0(0.),weight(0.),continuous_sfr(true){
			if(RSigmaZ==0.)RSigmaZ=RSigma;
			if(RD_IO  ==0.)RD_IO=RD;
			if(sigmaR0==0.)sigmaR0=sigmaR;
			if(sigmaZ0==0.)sigmaZ0=sigmaZ;
			if(weight==0.) continuous_sfr=true;
			else continuous_sfr=false;
			Norm = .5/(8.*PI*PI*PI);
		}
	void set_local_norm(bool r){local_norm=r;}
	void use_original_radius(bool r){use_original_radius_for_drift=r;}
	void set_average_radius(bool r){use_average_radius=r;}
	void set_sfr(SFR_base *sfr1){sfr=sfr1;}
	void set_mr(metal_relation *mr1){mr=mr1;}
	void setParams(VecDoub params){
		RD = params[0];		RSigma = params[1];
		sigmaR=params[2];		sigmaZ=params[3];
		if(params.size()>4){
			tau_T=params[4]; tau_m=params[5];
			RSigmaZ=params[6]; RD_IO=params[7];
			sigmaR0=params[8]; sigmaZ0=params[9];
			weight=params[10];
		}
		if(RSigmaZ==0.)RSigmaZ=RSigma;
		if(RD_IO  ==0.)RD_IO=RD;
		if(sigmaR0==0.)sigmaR0=sigmaR;
		if(sigmaZ0==0.)sigmaZ0=sigmaZ;
		if(weight==0.) continuous_sfr=true;
		else continuous_sfr=false;
	}
	double Rd(double tau); // scalelength for inside-out formation
	double SigmaRthick(double Rc, double age);
	double SigmaZthick(double Rc, double age);
	double fullDF_actions(const VecDoub& Actions, double age, double RcP);
	double chemDF_actions(const VecDoub& Actions, double F);
	double chemDF_real(const VecDoub& X, double F);
	double full_DF(const VecDoub& X, double age, double RcP);
	double full_DF_Z(const VecDoub& X, double age, double Z);
	void setStandardSolar(VecDoub SS){StandardSolar=SS;}
};
#endif
