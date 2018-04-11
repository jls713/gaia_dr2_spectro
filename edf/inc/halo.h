#ifndef HALOH
#define HALOH
// ===========================================================================
// Halo EDF
// ===========================================================================
#include "edf_base.h"
// ===========================================================================

class halo_edf: public edf{
protected:
	const double HALOWEIGHT = 4.27e-04*0.158655*0.25;
	double J0 = 180., phicoeff=0.68, zcoeff=.7, halo_slope=3.;
public:
	double HW, FHALO, SIGFHALO, min_halo_age=12.;
	halo_edf(Potential_JS *Pot, Actions_AxisymmetricFudge_InterpTables *Act,VecDoub params):edf(Pot,Act),HW(params[0]),FHALO(params[1]),SIGFHALO(params[2]){}
	void setParams(VecDoub params){
		HW = params[0];		FHALO = params[1];
		SIGFHALO=params[2];
		if(params.size()>3){
			if(params[3]!=0.) J0=params[3];
			if(params[4]!=0.) phicoeff=params[4];
			if(params[5]!=0.) zcoeff=params[5];
			if(params[6]!=0.) halo_slope=params[6];
			if(params[7]!=0.) min_halo_age=params[7];
		}
	}
	double chemDF_actions(const VecDoub& Actions, double F);
	double chemDF_real(const VecDoub& X, double F);
	double halo_metal(double F);
	double halo_actions(const VecDoub& Actions);
	double haloDF_nometal_real(const VecDoub& X);
};
// ===========================================================================
#endif
