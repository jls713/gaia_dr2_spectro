#ifndef THIN_DISCH
#define THIN_DISCH
// ===========================================================================
// Thin EDF
// ===========================================================================
#include "edf_base.h"
#include "sfr.h"
#include "metal_relation.h"
// ===========================================================================

class thin_disc_edf: public edf{
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
	bool local_norm=false; // if true, the SFR is the surface density SFR at the Sun (averaged over the column), else it is the global star formation rate per unit mass
public:

	double RD,RSigma,sigmaR,sigmaZ,tau_T,tau_1,BETA_R,BETA_Z;
	double RD_IO=0.,RSigmaZ=0.,tau_1_z=0.,tau_sigma=0.;

	thin_disc_edf(Potential_JS *Pot, Actions_AxisymmetricFudge_InterpTables *Act,metal_relation *mr,SFR_base *sfr,RadialMigration *rm, VecDoub params,int order=10)
		:edf(Pot,Act),mr(mr),sfr(sfr),rm(rm),GLI(order)
		,RD(params[0]),RSigma(params[1]),sigmaR(params[2])
		,sigmaZ(params[3]),tau_T(params[4]),tau_1(params[5])
		,BETA_R(params[6]),BETA_Z(params[7]),RD_IO(0.),RSigmaZ(0.),tau_1_z(0.),tau_sigma(params[4]){
			if(RD_IO==0.)  RD_IO=RD;
			if(RSigmaZ==0.)RSigmaZ=RSigma;
			if(tau_1_z==0.)tau_1_z=tau_1;
			Norm = .5/(8.*PI*PI*PI);
		}
	void set_local_norm(bool r){local_norm=r;}
	void set_average_radius(bool r){use_average_radius=r;}
	void set_sfr(SFR_base *sfr1){sfr=sfr1;}
	void set_mr(metal_relation *mr1){mr=mr1;}

	void setParams(VecDoub params){
		RD = params[0];		RSigma = params[1];
		sigmaR=params[2];		sigmaZ=params[3];
		if(params.size()>4){
			tau_T=params[4];		tau_1=params[5];
			BETA_R=params[6];		BETA_Z=params[7];
			RD_IO=params[8]; 		RSigmaZ=params[9];
			tau_1_z=params[10];     tau_sigma=params[11];
		}
		if(RD_IO==0.)  RD_IO=RD;
		if(RSigmaZ==0.)RSigmaZ=RSigma;
		if(tau_1_z==0.)tau_1_z=tau_1;
		if(tau_sigma==0.)tau_sigma=tau_T;
	}
	double Rd(double tau); // scalelength for inside-out formation
	double SigmaRthin(double Rc, double age){ // sigmaR is age dep. for thin disc
		return sigmaR*pow((age+tau_1)/(tau_sigma+tau_1),BETA_R)*exp((StandardSolar[0]-Rc)/RSigma);
	}
	double SigmaZthin(double Rc, double age){ // sigmaZ is age dep. for thin disc
		return sigmaZ*pow((age+tau_1_z)/(tau_sigma+tau_1_z),BETA_Z)*exp((StandardSolar[0]-Rc)/RSigmaZ);
	}
	double fullDF_actions(const VecDoub& Actions, double age, double RcP);
	double chemDF_actions(const VecDoub& Actions, double F);
	double chemDF_real(const VecDoub& X, double F);
	double full_DF(const VecDoub& X, double age, double RcP);
	double full_DF_Z(const VecDoub& X, double age, double Z);
	void setStandardSolar(VecDoub SS){StandardSolar=SS;}
};

///===========================================================================

#endif
