#include "sfr.h"

static double SFRGlobal(double tau, void *p){
	SFR_st *P = (SFR_st *) p;
	return exp(tau/P->T0-P->tau_turn/(P->tau_m-tau));
}

SFR::SFR(double T0, double tau_turn, double tau_m, double totalMass):SFR_base(tau_m,tau_turn),T0(T0),totalMass(totalMass){
	SFR_st P(T0,tau_m,tau_turn);
	SFRNorm = GaussLegendreQuad(&SFRGlobal,0.,tau_m-0.001,&P)/totalMass;
}

void SFR::reset(VecDoub params){
	T0=params[0]; tau_turn=params[1]; tau_m=params[2]; totalMass = params[3];
	SFR_st P(T0,tau_m,tau_turn);
	SFRNorm = GaussLegendreQuad(&SFRGlobal,0.,tau_m-0.001,&P)/totalMass;
}

void SFR_interp::reset(VecDoub vals1, double taum){
	tau_m=tau_m;vals=vals1;
	double delta_x = tau_m/(n_nodes_i-1.);
	double norm = 0.;
	for(int i=0;i<n_nodes_i;++i)
		norm+=vals[i]*delta_x;
	for(int i=0;i<n_nodes_i;++i)
		vals[i]/=norm;
	delete Int; Int = new interpolator(nodes,vals);
}
