#ifndef METAL_RELATIONH
#define METAL_RELATIONH
// ===========================================================================
// Metallicity with radius and time relations
// ===========================================================================
#include "utils.h"
#include "coordtransforms.h"
#include "tables_aa.h"
#include "rad_mig.h"

// Base class
class metal_relation{
protected:
public:
	double FEHStart, FEHGrad, FEHOffset, TAUF, tau_m;
	metal_relation(double FEHStart, double FEHGrad, double FEHOffset, double TAUF, double tau_m): FEHStart(FEHStart), FEHGrad(FEHGrad), FEHOffset(FEHOffset), TAUF(TAUF), tau_m(tau_m){}
	void reset(VecDoub params){
		FEHStart=params[0]; FEHGrad=params[1]; FEHOffset=params[2];
		TAUF=params[3]; tau_m=params[4];
	}
	inline double FStart(void){return FEHStart;}
	virtual double FeH(double age, double Rc){return 0.;}
	virtual double AgeFromMetalRadius(double Z, double R){return 0.;}
	virtual double RadiusFromMetal(double age, double FeH){return 0.;}
	virtual double dZdRc(double Rc, double tau){return 0.;}
	virtual double FR(double RcP){return 0.;}
	inline double MaxAge(double Z){return AgeFromMetalRadius(Z,0.);}
	double tauDF(double RcPrimed, double F, double *AGE);
	double MetalDF(double Lz, double FeH, double age, double *RCPrimed, double *LzPrimed, Actions_AxisymmetricFudge_InterpTables *Act, RadialMigration *rm,double RD,double Jz=0.);
};
// ===========================================================================

// Tanh relation where FeH(R)= F Tanh (dFdR*(R-R(0))/FeH_initial)
class tanh_metal_relation: public metal_relation{
private:
public:
	tanh_metal_relation(double FEHStart, double FEHGrad, double FEHOffset, double TAUF, double tau_m):metal_relation(FEHStart,FEHGrad,FEHOffset,TAUF,tau_m){};
	inline double FeH(double age, double Rc){
		double FeHL0 = FEHStart*tanh(FEHGrad*(Rc-FEHOffset)/FEHStart);
		return FeHL0+(FeHL0-FEHStart)*(tanh((tau_m-age)/TAUF)-1);
	}
	inline double AgeFromMetalRadius(double Z, double R){
		double FR0=FEHStart*tanh(FEHGrad*(R-FEHOffset)/FEHStart);
		return tau_m-TAUF*atanh((Z-FEHStart)/(FR0-FEHStart));
	}
	inline double RadiusFromMetal(double age, double FeH){
		// Will return less than zero if Z>Z_max(age)
		double g = (tanh((tau_m-age)/TAUF)-1);
		double fL0=(FeH/FEHStart+g)/(1+g);
		if(fabs(fL0)>1.) return -10.;
		return (FEHStart/FEHGrad)*atanh(fL0)+FEHOffset;
	}
	inline double dZdRc(double Rc, double tau){
		double g = tanh((tau_m-tau)/TAUF);
		double cc = cosh(FEHGrad*(Rc-FEHOffset)/FEHStart);
		return FEHGrad*g/cc/cc;
	}
	inline double FR(double Rc){
		return FEHStart*tanh(FEHGrad*(Rc-FEHOffset)/FEHStart);
	}
};
// ===========================================================================

// Exp relation where FeH(R)= F (1-Exp(-dFdR*(R-R(0))/FeH_initial)
class exp_metal_relation: public metal_relation{
private:
public:
	exp_metal_relation(double FEHStart, double FEHGrad, double FEHOffset, double TAUF, double tau_m):metal_relation(FEHStart,FEHGrad,FEHOffset,TAUF,tau_m){};
	inline double FeH(double age, double Rc){
		double FeHL0 = FEHStart*(1-exp(-FEHGrad*(Rc-FEHOffset)/FEHStart));
		return FeHL0+(FeHL0-FEHStart)*(tanh((tau_m-age)/TAUF)-1);
	}
	inline double AgeFromMetalRadius(double Z, double R){
		double FR0 = FEHStart*(1-exp(-FEHGrad*(R-FEHOffset)/FEHStart));
		return tau_m-TAUF*atanh((Z-FEHStart)/(FR0-FEHStart));
	}
	inline double RadiusFromMetal(double age, double FeH){
		// Will return less than zero if Z>Z_max(age)
		double g = (tanh((tau_m-age)/TAUF)-1);
		double fL0=(FeH/FEHStart+g)/(1+g);
		if(fL0>1.) return -10.;
		return -(FEHStart/FEHGrad)*log(1.-fL0)+FEHOffset;
	}
	inline double dZdRc(double Rc, double tau){
		double g = tanh((tau_m-tau)/TAUF);
		return FEHGrad*g*exp(-FEHGrad*(Rc-FEHOffset)/FEHStart);
	}
	inline double FR(double Rc){
		return FEHStart*(1-exp(-FEHGrad*(Rc-FEHOffset)/FEHStart));
	}

};
// Linear type relation where FeH(R,t)=log(exp(FEH_i)-(age-tau_m)/tau_m*(dFdR*(R-R(0)-exp(FEH_i))
// TAUF doesn't do anything
class linear_metal_relation: public metal_relation{
private:
public:
	linear_metal_relation(double FEHStart, double FEHGrad, double FEHOffset, double TAUF, double tau_m):metal_relation(exp(FEHStart),FEHGrad,FEHOffset,0.,tau_m){};
	inline double FeH(double age, double Rc){
		double FeHL0 = FEHGrad*(Rc-FEHOffset);
		if(FeHL0<FEHStart)
			return log(FEHStart);
		return log(FEHStart+(age-tau_m)*(FeHL0-FEHStart)/tau_m);
	}
	inline double AgeFromMetalRadius(double Z, double R){
		double FR0 = FEHGrad*(R-FEHOffset)-FEHStart;
		return ((exp(Z)-FEHStart)/FR0+1.)*tau_m;
	}
	inline double RadiusFromMetal(double age, double FeH){
		// Will return less than zero if Z>Z_max(age)
		double g = (age/tau_m-1.);
		double f = ((exp(FeH)-FEHStart)/g+FEHStart)/FEHGrad+FEHOffset;
		if(f<0.) return -10.;
		return f;
	}
	inline double dZdRc(double Rc, double tau){
		double FeHL0 = FEHGrad*(Rc-FEHOffset);
		double g = (tau/tau_m-1.)*FEHGrad;
		return g/(FEHStart+(tau/tau_m-1.)*(FeHL0-FEHStart));
	}
	inline double FR(double Rc){
		double FeHL0 = FEHGrad*(Rc-FEHOffset);
		return log(FeHL0);
	}
};
// ===========================================================================

// Simple Alpha/Fe as a function of radius and time relation
class alpha_relation{
protected:
public:
	double a_init, a_finals, a_grad, ts0, tsup, Rs, a0, a1, ta;
	alpha_relation(double a_init,double a_finals, double a_grad, double ts0, double tsup, double Rs)
		: a_init(a_init), a_finals(a_finals), a_grad(a_grad),ts0(ts0), Rs(Rs) {
	}
	void reset(VecDoub params){
		a_init=params[0]; a_finals=params[1]; a_grad=params[2];
		ts0=params[3]; tsup=params[4]; Rs=params[5];
	}
	inline double AStart(void){return a_init;}
	double Alpha(double age, double Rc){
    		double ts = ts0*exp(Rc/Rs);
    		double ta = 12.-tsup;
    		double a_f=a_grad*(Rc-conv::StandardSolar[0])+a_finals;
   	 	double a1 = (a_init-a_f)/(tanh(ta/ts)+tanh(tsup/ts));
    		double a0 = a_init-a1*tanh(tsup/ts);
    		return a0+a1*tanh((age-ta)/ts);
	}
	double AgeFromAlphaRadius(double A, double R){return 0.;}
	double RadiusFromAlpha(double age, double A){return 0.;}
};
// ===========================================================================
#endif
