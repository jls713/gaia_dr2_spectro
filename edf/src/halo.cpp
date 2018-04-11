/*=============================================================================
  Halo DF
=============================================================================*/
#include "halo.h"
//=============================================================================

double halo_edf::halo_metal(double F){
	return exp(-(F-FHALO)*(F-FHALO)/2./SIGFHALO/SIGFHALO)/sqrt(2.*PI*SIGFHALO*SIGFHALO);
}

double halo_edf::halo_actions(const VecDoub& Actions){
	double JR = Actions[0],Jp = Actions[1],Jz = Actions[2];
	return pow(J0+JR+phicoeff*fabs(Jp)+zcoeff*Jz,-halo_slope);
}

double halo_edf::haloDF_nometal_real(const VecDoub& X){
	if(pot->H(X)>0.) return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	return HALOWEIGHT*HW*halo_actions(actions);
}

double halo_edf::chemDF_actions(const VecDoub& A, double F){
	return HALOWEIGHT*HW*halo_actions(A)*halo_metal(F);
}

double halo_edf::chemDF_real(const VecDoub& X, double F){
	if(pot->H(X)>0.) return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	double P = chemDF_actions(actions,F);
	if(P<0.)std::cerr<<"DF negative: Z="<<F<<std::endl;
	if(P!=P)std::cerr<<"DF nan: Z="<<F<<std::endl;
	return P;
}
