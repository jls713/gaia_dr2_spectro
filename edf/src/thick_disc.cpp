/*=============================================================================
  Thick Disc DF
=============================================================================*/
#include "thick_disc.h"
//=============================================================================

double thick_disc_edf::Rd(double tau){
	return RD+(RD_IO-RD)*(tau-tau_T)/(tau_m-tau_T);
}

double thick_disc_edf::SigmaRthick(double Rc, double age){
		return (sigmaR+(sigmaR0-sigmaR)*(age-tau_T)/(tau_m-tau_T))*exp((StandardSolar[0]-Rc)/RSigma);
	}
double thick_disc_edf::SigmaZthick(double Rc, double age){
	return (sigmaZ+(sigmaZ0-sigmaZ)*(age-tau_T)/(tau_m-tau_T))*exp((StandardSolar[0]-Rc)/RSigmaZ);
}

double thick_disc_edf::fullDF_actions(const VecDoub& Actions, double age, double RcP){
	double Rc = Actions[3], JR = Actions[0], Jz = Actions[2];
	double Lz = Actions[1], kappa = Actions[4], nuz = Actions[5];
	if(use_average_radius){
		Rc=(Rc+RcP)/2.;
	}
	if(usecombinationLzJzforSigma){
		// find Rc for Lz & Jz combo
		Rc = pot->R_L(Lz+Jz,Rc);
	}
	double sigmar=SigmaRthick(Rc,age);
	double sigmar2=pow(sigmar,2);
	double sigmaz=SigmaZthick(Rc,age);
	double sigmaz2=pow(sigmaz,2);
	double R_d = Rd(age);
	double Sd = Norm*(1+tanh(Lz/L0))*exp(-RcP/R_d);
	if(!local_norm) Sd/=(R_d*R_d);
	else Sd*=exp(StandardSolar[0]/R_d);
	if(continuous_sfr) Sd*=sfr->SFR_w(age);
	else Sd*=weight;
	return exp(-kappa*JR/sigmar2-nuz*Jz/sigmaz2)*kappa*nuz/(sigmar2*sigmaz2)*Sd;
}

double thick_disc_edf::chemDF_actions(const VecDoub& Actions, double FeH){
	// for the thick disc we use the delta function to perform the age
	// integral, and then integrate over birth ang. mom. using Gaussian quadrature.
	// This is because dF/dL -> 0 at t=tau_m

	double Lz=Actions[1], Jz=Actions[2], Rc=Actions[3];

	double dR, Rplus, Rminus, tau, dFdt, SigF, Partial, OmkR, Lp, R_d, DV;
	double fplus=0., fminus=0., ss=0.;
	// We find the maximum radius at some metallicity by finding the radius at
	// F = FeH, tau = 10
	double RCmax = mr->RadiusFromMetal(tau_T,FeH);
	if(RCmax<0.) return 0.;
	double SIGTHIN = rm->SigmaF(tau_m,Rc,Jz);
	double b=Rc+3.*SIGTHIN/220., a = Rc-3.*SIGTHIN/220.;
	if(a<0.)a=0.001;
	if(b>RCmax)b=RCmax;
	double Rm = (b+a)/2., Rr = (b-a)/2.;
	ss=0.;
	int order=GLI.get_order()/2; VecDoub abweight;
	for(int i=0;i<order;i++){
		abweight=GLI.abscissa_and_weight_half(i);
		dR=Rr*abweight[0];
		Rminus=Rm-dR; Rplus=Rm+dR;
		if(Rminus>0.){
		Lp =ActionCalculator->L_circ_current(Rminus);
		dFdt=mr->tauDF(Rminus,FeH,&tau);
		SigF=rm->SigmaF(tau,Rminus,Jz);
		if(tau>tau_T and tau<tau_m and dFdt>0.){ //thick disc
			R_d = Rd(tau);
			DV = rm->DriftVelocity(tau,Rminus,Jz,
			            use_original_radius_for_drift?R_d:0.);
			Partial=fullDF_actions(Actions,tau,Rminus);
			OmkR=Rminus/rm->LzPrimeNormalization(Lp+DV, tau, Jz);
			fminus=OmkR*Partial*GFunction(Lz-Lp-DV, SigF)/dFdt;
			if(fminus!=fminus) fminus=0.;
		}
		}
		else fminus=0.;

		Lp = ActionCalculator->L_circ_current(Rplus);
		dFdt=mr->tauDF(Rplus,FeH,&tau);
		SigF=rm->SigmaF(tau,Rplus,Jz);
		if(tau>tau_T and tau<tau_m and dFdt>0.){ // thick disc
			R_d = Rd(tau);
			DV = rm->DriftVelocity(tau,Rplus,Jz,
			            use_original_radius_for_drift?R_d:0.);
			Partial=fullDF_actions(Actions,tau,Rplus);
			OmkR=Rplus/rm->LzPrimeNormalization(Lp+DV, tau, Jz);
			fplus=OmkR*Partial*GFunction(Lz-Lp-DV, SigF)/dFdt;
			if(fplus!=fplus) fplus=0.;
		}
		else fplus=0.;

		ss=ss+abweight[1]*(fplus+fminus);
	}
	return ss*Rr;
}

double thick_disc_edf::chemDF_real(const VecDoub& X, double F){
	if(F<mr->FStart())return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	double P = chemDF_actions(actions,F);
	if(P<0.)std::cerr<<"DF negative: Z="<<F<<std::endl;
	if(P!=P)std::cerr<<"DF nan: Z="<<F<<std::endl;
	return P;
}
double thick_disc_edf::full_DF(const VecDoub& X, double age, double RcP){
	double FeH = mr->FeH(age,RcP);
	if(FeH<mr->FStart())return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	double P = fullDF_actions(actions,age,RcP);
	if(P<0.)std::cerr<<"DF negative: Z="<<FeH<<std::endl;
	if(P!=P)std::cerr<<"DF nan: Z="<<FeH<<std::endl;
	return P;
}
double thick_disc_edf::full_DF_Z(const VecDoub& X, double age, double FeH){
	if(FeH<mr->FStart() or age<tau_T)return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	double RcP = mr->RadiusFromMetal(age,FeH);
	double P = fullDF_actions(actions,age,RcP);
	double dRdz = fabs(1./mr->dZdRc(RcP, age));
	if(P<0.)std::cerr<<"DF negative: Z="<<FeH<<std::endl;
	if(P!=P)std::cerr<<"DF nan: Z="<<FeH<<std::endl;
	return P*dRdz;
}
