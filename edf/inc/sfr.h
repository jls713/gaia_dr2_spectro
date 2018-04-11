#ifndef SFRH
#define SFRH
// ===========================================================================
// Star Formation Rate
// ===========================================================================
#include "utils.h"
// #include "GSLInterface.h"
// ===========================================================================

// Used to integrate SFR(tau)
struct SFR_st{
	double T0, tau_m, tau_turn, totalMass;
	SFR_st(double T0, double tau_m, double tau_turn)
		:T0(T0), tau_m(tau_m), tau_turn(tau_turn){}
};


class SFR_base{
public:
	double tau_m, tau_turn;
	SFR_base(double taum, double tau_turn):tau_m(taum), tau_turn(tau_turn){}
	virtual double SFR_w(double tau)=0;
	virtual void reset(VecDoub params)=0;
	virtual void reset_mass(double mass)=0;
};

class SFR: public SFR_base{
public:
	double T0, SFRNorm, totalMass;
	SFR(double T0, double tau_turn, double tau_m, double totalMass);
	double SFR_w(double tau){
		return exp(tau/T0-tau_turn/(tau_m-tau))/SFRNorm;
	}
	void reset(VecDoub params);
	void reset_mass(double m){
		std::cout<<"Resetting disc mass, SFRNorm=";
		SFRNorm*=totalMass/m;totalMass=m;
		std::cout<<SFRNorm<<std::endl;
	}
};

class SFR_interp: public SFR_base{
public:
	int n_nodes_i;
	VecDoub nodes, vals;
	interpolator *Int;
	SFR_interp(VecDoub vals1, int n_nodes,double tau_m)
		:SFR_base(tau_m,0.),n_nodes_i(n_nodes),nodes(VecDoub(n_nodes,0.)),vals(vals1){
			double delta_x = tau_m/(n_nodes-1.);
			for(int i=0;i<n_nodes;++i)
				nodes[i]=i*delta_x;
			double norm = 0.;
			for(int i=0;i<n_nodes;++i)
				norm+=vals[i]*delta_x;
			for(int i=0;i<n_nodes;++i)
				vals[i]/=norm;
			Int=new interpolator(nodes,vals);
		}
	double SFR_w(double tau){
		return Int->interpolate(tau);
	}
	void reset(VecDoub vals){
		return reset(vals,tau_m);
	}
	void reset(VecDoub vals, double tau_m);
	void reset_mass(double m){return;}
};
#endif
