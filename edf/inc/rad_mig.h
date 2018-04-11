#ifndef RMH
#define RMH
// ===========================================================================
// Radial Migration Kernel
// ===========================================================================

class RadialMigration{
private:
public:
	double SIGTHIN, tau_m, Vcflat, GAMMAT, Rd;
	double R0;
	double coeff;
	double Jz0;
	double sigR,inv2sigR;
	RadialMigration(double SIGTHIN, double tau_m, double GAMMAT, double Rd, Potential_JS *Pot,double R0, double Jz0=0.,double sigR=0.):SIGTHIN(SIGTHIN),tau_m(tau_m),Vcflat(Pot->Vc(R0)),GAMMAT(GAMMAT),Rd(Rd),R0(R0),Jz0(Jz0),sigR(sigR){
		coeff = -pow(SIGTHIN,2.)*0.5/tau_m/Rd/Vcflat;
		inv2sigR=1./sigR/sigR;
	}
	void reset(VecDoub params,Potential_JS *Pot,double R0, double Jz04=0.,double sigR4=0.){
		SIGTHIN = params[0];	tau_m = params[1];
		GAMMAT=params[2];		Rd=params[3];
		if(params.size()>4) Jz0=params[4];
		Vcflat = Pot->Vc(R0);
		coeff = -pow(SIGTHIN,2.)*0.5/tau_m/Rd/Vcflat;
		sigR=sigR4;inv2sigR=1./sigR/sigR;
		Jz0=Jz04;
	}
	inline double SigmaF(double age,double Rc, double Jz=0.){
		if(age==0.)age+=1e-6;
		double fac=1.;
		if(sigR)
			fac = sqrt(Rc/R0*exp(-(Rc-R0)*(Rc-R0)*.5*inv2sigR));
		return (Jz0>0.?exp(-Jz/Jz0):1.)*fac*SIGTHIN*pow((age/tau_m),GAMMAT);
	}
	inline double DriftVelocity(double age,double Rc, double Jz=0., double RD=0.){
		double fac=1.;
		if(sigR)
			fac *= Rc/R0*exp(-(Rc-R0)*(Rc-R0)*.5*inv2sigR);
		if(RD)
			fac *= Rd/RD;
		return (Jz0>0.?exp(-2.*Jz/Jz0):1.)*fac*coeff*age;
	}
	inline double LzPrimeNormalization(double LzP, double age, double Jz=0.){
		double P = 0.5*(1+erf(LzP/sqrt(2.)/SigmaF(age, LzP/Vcflat ,Jz)));
		if(P>0.) return P;
		else return 1e-20;
	}
};
// ===========================================================================
#endif
