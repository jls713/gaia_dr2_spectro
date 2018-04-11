#include "metal_relation.h"

double metal_relation::MetalDF(double Lz, double FeH, double age, double *RCPrimed, double *LzPrimed, Actions_AxisymmetricFudge_InterpTables *Act, RadialMigration *rm, double RD,double Jz){

	double g = (tanh((tau_m-age)/TAUF)-1);
	double fL0=(FeH/FEHStart+g)/(1+g);
	if(fabs(fL0)>1.) return 0.;

	double RC0 = RadiusFromMetal(age,FeH);
	if(RC0<0.) return 0.;

	double LZ0 = Act->L_circ_current(RC0);
	(*RCPrimed)=RC0;
	(*LzPrimed)=LZ0;
	double SigF=rm->SigmaF(age,RC0,Jz);

	double kappa_h, Omega_h;
	// Pot->getfreqs(RC0,&kappa_h,&nu_h,&Omega_h);
	VecDoub freqs = Act->PotentialFreq(RC0);
	kappa_h = freqs[0]; Omega_h = freqs[2];
	double dLdRc=kappa_h*kappa_h*RC0/2./Omega_h;
	double dFdL=fabs(dZdRc(RC0,age)/dLdRc);
	double DV = rm->DriftVelocity(age,RC0,Jz);
	return GFunction(Lz-LZ0-DV,SigF)/dFdL;

}
double metal_relation::tauDF(double RcPrimed, double F, double *AGE){
	double fr=FR(RcPrimed);
	double g = (F-fr)/(fr-FEHStart)+1.;
	if(fabs(g)>1.){ (*AGE)=1e99;return -10.;}
	double tau = tau_m-TAUF*atanh(g);
	(*AGE)=tau;
	double S = cosh((tau_m-tau)/TAUF);
	return ((fr-FEHStart)/(TAUF*S*S));
}
