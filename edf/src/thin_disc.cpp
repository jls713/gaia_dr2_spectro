#include "thin_disc.h"
/*=============================================================================
 	DISC DISTRIBUTION FUNCTION

	The disc DF is composed of a thin and thick disc. The dynamical part consists of a
	superposition of QI for the thin disc and a single QI for the thick disc.
	The metallicity part consists of multiplying the DF by a Gaussian in angular momentum
	to represent the churning with the mean given by the birth angular momentum which is
	a function of the metallicity of the star, and its age.

=============================================================================*/

double thin_disc_edf::Rd(double tau){
	return RD+(RD_IO-RD)*tau/tau_T;
}

double thin_disc_edf::fullDF_actions(const VecDoub& Actions, double age, double RcP){
	double Rc = Actions[3], JR = Actions[0], Jz = Actions[2];
	double Lz = Actions[1], kappa = Actions[4], nuz = Actions[5];
	if(use_average_radius){
		Rc=(Rc+RcP)/2.;
	}
	if(usecombinationLzJzforSigma){
		// find Rc for Lz & Jz combo
		Rc = pot->R_L(Lz+Jz,Rc);
	}
	double sigmar=SigmaRthin(Rc,age);
	double sigmar2=pow(sigmar,2);
	double sigmaz=SigmaZthin(Rc,age);
	double sigmaz2=pow(sigmaz,2);
	double R_d = Rd(age);
	double Sd = Norm*(1+tanh(Lz/L0))*exp(-RcP/R_d);
	if(!local_norm) Sd/=(R_d*R_d);
	else Sd*=exp(StandardSolar[0]/R_d);
	// double Sd = 0.5*(1+tanh(Lz/L0))*exp(-RcP/R_d)/(8.*PI*PI*PI*R_d*R_d);
	return sfr->SFR_w(age)*exp(-kappa*JR/sigmar2-nuz*Jz/sigmaz2)*kappa*nuz/(sigmar2*sigmaz2)*Sd;
}

double thin_disc_edf::chemDF_actions(const VecDoub& Actions, double FeH){
	// for the thin disc we use the delta function to perform the birth ang. mom.
	// integral, and then integrate over tau using Gaussian quadrature.
	// This is because dF/dt -> 0 at t=0
	double Lz=Actions[1],Jz=Actions[2];
	// Thin
	double tauMIN=0.;double tauMAX=mr->tau_m;
	double agem,ager;
	double da,ageplus,ageminus,ss=0.0, RcPrimed=0, LzPrimed=0, kappaP, OmegaP, fplus=0., fminus=0., OmkR, MDF, R_d, DV;
	double a, b; VecDoub ff;
	a=tauMIN;
	if(tauMAX<tau_T) b=tauMAX;
	else b=tau_T;
	agem=0.5*(b+a);
	ager=0.5*(b-a);
	ss=0.0;
	int order=GLI.get_order()/2; VecDoub abweight;
	for(int i=0;i<order;i++){
		abweight=GLI.abscissa_and_weight_half(i);
		da=ager*abweight[0];
		ageplus=agem+da; ageminus=agem-da;
		RcPrimed=0;LzPrimed=0;fplus=0.;fminus=0.;
		R_d = Rd(ageplus);
		MDF=mr->MetalDF(Lz, FeH, ageplus, &RcPrimed, &LzPrimed, ActionCalculator, rm, R_d,Jz);
		DV = rm->DriftVelocity(ageplus,RcPrimed,Jz,R_d);
		if(MDF>0.){
			ff = ActionCalculator->PotentialFreq(RcPrimed);
			kappaP = ff[0];OmegaP = ff[2];
			OmkR=2.*OmegaP/kappaP/kappaP/rm->LzPrimeNormalization(LzPrimed+DV, ageplus,Jz);
			fplus=OmkR*fullDF_actions(Actions, ageplus, RcPrimed)*MDF;
			if(fplus!=fplus) fplus=0.;
		}
		R_d = Rd(ageminus);
		MDF=mr->MetalDF(Lz, FeH, ageminus, &RcPrimed, &LzPrimed,ActionCalculator,rm, R_d,Jz);
		DV = rm->DriftVelocity(ageminus,RcPrimed,Jz,R_d);
		if(MDF>0.){
			ff = ActionCalculator->PotentialFreq(RcPrimed);
			kappaP = ff[0];OmegaP = ff[2];
			OmkR=2.*OmegaP/kappaP/kappaP/rm->LzPrimeNormalization(LzPrimed+DV, ageminus,Jz);
			fminus=OmkR*fullDF_actions(Actions, ageminus, RcPrimed)*MDF;
			if(fminus!=fminus) fminus=0.;
		}
		ss=ss+abweight[1]*(fplus+fminus);
	}
	return ss*ager;
}

double thin_disc_edf::chemDF_real(const VecDoub& X, double F){
	if(F<mr->FStart())return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	double P = chemDF_actions(actions,F);
	if(P<0.)std::cerr<<"DF negative: Z="<<F<<std::endl;
	if(P!=P)std::cerr<<"DF nan: Z="<<F<<std::endl;
	return P;
}
double thin_disc_edf::full_DF(const VecDoub& X, double age, double RcP){
	double FeH = mr->FeH(age,RcP);
	if(FeH<mr->FStart() or age>tau_T)return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	double P = fullDF_actions(actions,age,RcP);
	if(P<0.)std::cerr<<"DF negative: Z="<<FeH<<std::endl;
	if(P!=P)std::cerr<<"DF nan: Z="<<FeH<<std::endl;
	return P;
}
double thin_disc_edf::full_DF_Z(const VecDoub& X, double age, double FeH){
	if(FeH<mr->FStart() or age>tau_T)return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	double RcP = mr->RadiusFromMetal(age,FeH);
	double P = fullDF_actions(actions,age,RcP);
	double dRdz = fabs(1./mr->dZdRc(RcP, age));
	if(P<0.)std::cerr<<"DF negative: Z="<<FeH<<std::endl;
	if(P!=P)std::cerr<<"DF nan: Z="<<FeH<<std::endl;
	return P*dRdz;
}
