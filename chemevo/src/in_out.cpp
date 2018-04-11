#include "in_out.h"
//=============================================================================
DoubleInfallInflow::DoubleInfallInflow(ModelParameters M,double present_rSFR){
	auto F = M.parameters["flows"]["inflow"];
    std::vector<std::string> vars = {"FastTimeScale","SlowTimeScale","Weight"};
    for(auto v:vars)
        if (F.find(v) == F.end()) {
            LOG(INFO)<<v<<" not found in parameters file\n";
            throw std::invalid_argument(v+" not found in parameters file");
        }
    F = M.parameters["fundamentals"];
    vars = {"PresentInfallRate","GasScaleLength"};
    for(auto v:vars)
        if (F.find(v) == F.end()) {
            LOG(INFO)<<v<<" not found in parameters file\n";
            throw std::invalid_argument(v+" not found in parameters file");
        }
	tf = M.parameters["flows"]["inflow"]["FastTimeScale"];
	ts = M.parameters["flows"]["inflow"]["SlowTimeScale"];
	weight = M.parameters["flows"]["inflow"]["Weight"];
	Rd = M.parameters["fundamentals"]["GasScaleLength"];
	pif=1.;
	double pif0 = (*this)(M.parameters["fundamentals"]["SolarRadius"],
	                      M.parameters["fundamentals"]["GalaxyAge"]);
	pif = (double)M.parameters["fundamentals"]["PresentInfallRate"]/pif0;
}
double DoubleInfallInflow::operator()(double R, double t, Grid *rSFR){
	return (exp(-t/ts)*(1.-weight)+exp(-t/tf)*weight)*pif*exp(-R/Rd);
}
//=============================================================================
SimpleGalacticFountain::SimpleGalacticFountain(ModelParameters M){
	auto F = M.parameters["flows"]["outflow"];
    std::vector<std::string> vars = {"Feject_in","Feject_out","TransitionRadius"};
    for(auto v:vars)
        if (F.find(v) == F.end()) {
            LOG(INFO)<<v<<" not found in parameters file\n";
            throw std::invalid_argument(v+" not found in parameters file");
        }
	feject_in=M.parameters["flows"]["outflow"]["Feject_in"];
	feject_out=M.parameters["flows"]["outflow"]["Feject_out"];
	transR = M.parameters["flows"]["outflow"]["TransitionRadius"];
	A=.5*(feject_out+feject_in);
	B=.5*(feject_out-feject_in);
	dR = 0.2;
}
double SimpleGalacticFountain::operator()(double R, double t){
	return A+B*tanh((R-transR)/dR);
}
//=============================================================================
double RadialFlow::beta_g(double R, double Rdown, double Rup, double t, double dt, Grid *rSFR, int*err){
	auto rhd = (R+Rdown)*.5;
	auto vd=flow_rate(rhd,t,rSFR);
	if(fabs(vd)>(R-Rdown)/dt)
		*err=1;
		// throw std::runtime_error("Courant condition not satisfied: v="
		//                          +std::to_string(vd)
		//                          +", DeltaR="+std::to_string(R-Rdown)
		//                          +", DeltaT="
		//                          +std::to_string(dt)
		//                          +", DeltaR/DeltaT="
		//                          +std::to_string((R-Rdown)/dt)
		//                          +"\n");
	// if(vd>0.)
	// 	LOG(INFO)<<"Outflows (v="<<std::to_string(vd)<<") at R="
	// 	                         <<std::to_string(R)<<",t="
	// 	                         <<std::to_string(t)<<"\n";
	return -2.*vd*(Rdown+R)/(Rup-Rdown)/(R+.5*(Rdown+Rup));
}

double RadialFlow::gamma_g(double R, double Rdown, double Rup, double t, double dt, Grid *rSFR, int*err){
	auto rhu = (Rup+R)*.5;
	auto vu=flow_rate(rhu,t,rSFR);
	if(fabs(vu)>(Rup-R)/dt)
		*err=1;
		// throw std::runtime_error("Courant condition not satisfied: v="
		//                          +std::to_string(vu)
		//                          +", DeltaR="+std::to_string(Rup-R)
		//                          +", DeltaT="
		//                          +std::to_string(dt)
		//                          +", DeltaR/DeltaT="
		//                          +std::to_string((Rup-R)/dt)
		//                          +"\n");
	// if(vu>0.)
	// 	LOG(INFO)<<"Outflows (v="<<std::to_string(vu)<<") at R="
	// 	                         <<std::to_string(R)<<",t="
	// 	                         <<std::to_string(t)<<"\n";
	return -2.*vu*(Rup+R)/(Rup-Rdown)/(R+.5*(Rdown+Rup));
}

double RadialFlow::dMdt(double mass, double massup, double R, double Rdown, double Rup, double t, double dt, Grid *rSFR, int*err){
	return -beta_g(R,Rdown,Rup,t,dt,rSFR,err)*mass+gamma_g(R,Rdown,Rup,t,dt,rSFR,err)*massup;
}
//=============================================================================
LinearRadialFlow::LinearRadialFlow(ModelParameters M,double present_rSFR){
	auto F = M.parameters["flows"]["radialflow"];
	std::string vars = "Gradient";
	if (F.find(vars) == F.end()) {
        LOG(INFO)<<vars<<" not found in parameters file\n";
        throw std::invalid_argument(vars+" not found in parameters file");
    }
	VGrad = M.parameters["flows"]["radialflow"]["Gradient"];
}
//=============================================================================
// PezzulliInflowRadialFlow::PezzulliInflowRadialFlow(ModelParameters M,std::shared_ptr<StarFormationRate> sfr):SFR(sfr){
// 	auto F = M.parameters["flows"]["radialflow"];
// 	std::string vars = "PezzulliAlpha";
// 	if (F.find(vars) == F.end()) {
// 		LOG(INFO)<<vars<<" not found in parameters file\n";
// 		throw std::invalid_argument(vars+" not found in parameters file");
//     }
// 	Alpha = M.parameters["flows"]["radialflow"]["PezzulliAlpha"];
// 	KSN = M.parameters["fundamentals"]["Kennicutt-Schmidt_Coeff"];
// 	double PresentGasDensity = M.parameters["fundamentals"]["PresentGasDensitySun"];
// 	double PresentSFR = M.parameters["fundamentals"]["PresentSFR"];
// 	A = PresentSFR/pow(PresentGasDensity,KSN);
// }
// double PezzulliInflowRadialFlow::sigmagas(double R, double t){
// 	return pow((*SFR)(R,t)/A,1./KSN);
// }
// double PezzulliInflowRadialFlow::sigmaeffdot(double R, double t){
// 	double sfr = (*SFR)(R,t);
// 	double dsfr = ((*SFR)(R,t+0.005)-(*SFR)(R,t-0.005))/0.01;
// 	if(t<0.005) dsfr = ((*SFR)(R,t+0.01)-(*SFR)(R,t))/0.01;
// 	return sfr*(1.+pow(sfr,1./KSN-2.)*dsfr/pow(A,1./KSN)/KSN);
// }
// double mu_integrand(double R, void *p){
// 	mu_st *P=(mu_st *) p;
//     return pow(R,1.+1./P->P->alpha())*P->P->sigmaeffdot(R,P->t);
// }
// double PezzulliInflowRadialFlow::mu(double R, double t){
// 	// radial mass flow rate in units M_sun Gyr^-1
// 	if(R==0.) return 0.;
// 	GaussLegendreIntegrator GL(150);
//     mu_st Q = {this,t};
//     return -2.*PI*pow(R,-1./Alpha)*GL.integrate(&mu_integrand,0.,R,&Q);
// }
// double PezzulliInflowRadialFlow::acc_rate(double R, double t){
// 	// surface accretion rate in units M_sun kpc^-2 Gyr^-1
// 	return -mu(R,t)/(2.*PI*Alpha*R*R);
// }
// double PezzulliInflowRadialFlow::flow_rate(double R, double t, Grid *rSFR){
// 	// radial flow rate in units kpc Gyr^-1
// 	if(R==0.) return 0.;
// 	return mu(R,t)/sigmagas(R,t)/(2.*PI*R);
// }
// double PezzulliInflowRadialFlow::operator()(double R, double t, Grid *rSFR){
// 	return acc_rate(R,t);
// }
//=============================================================================
template<class T>
PezzulliInflowRadialFlow_rSFR<T>::PezzulliInflowRadialFlow_rSFR(ModelParameters M,double present_rSFR){
	auto F = M.parameters["flows"]["radialflow"];
	std::string vars = "PezzulliAlpha";
	if (F.find(vars) == F.end()) {
		LOG(INFO)<<vars<<" not found in parameters file\n";
		throw std::invalid_argument(vars+" not found in parameters file");
    }
	Alpha = M.parameters["flows"]["radialflow"]["PezzulliAlpha"];
	KSN = M.parameters["fundamentals"]["Kennicutt-Schmidt_Coeff"];
	double PresentGasDensity = M.parameters["fundamentals"]["PresentGasDensitySun"];
	A = present_rSFR/pow(PresentGasDensity,KSN);
}
double find_SFR(double R, double t, Grid* rSFR){
	// Extrapolate off edge of grid if necessary
	double sfr=0.;
	// Extrapolate off edge of grid
	if(R>rSFR->grid_radial().back()){
		// return 0.;
		sfr = rSFR->log_extrapolate_high(R,t);
		// auto NR = rSFR->grid_radial().size();
		// double rmax = rSFR->grid_radial()[NR-1];
		// double rmaxm = rSFR->grid_radial()[NR-2];
		// sfr = exp(log((*rSFR)(rmax,t))+(log((*rSFR)(rmax,t))-log((*rSFR)(rmaxm,t)))/(rmax-rmaxm)*(R-rmax));
	}
	else if(R<rSFR->grid_radial().front()){
		sfr = rSFR->log_extrapolate_low(R,t);
		// double rmax = rSFR->grid_radial()[1];
		// double rmaxm = rSFR->grid_radial()[0];
		// sfr = exp(log((*rSFR)(rmaxm,t))+(log((*rSFR)(rmax,t))-log((*rSFR)(rmaxm,t)))/(rmax-rmaxm)*(R-rmaxm));
	}
	else
		sfr = (*rSFR)(R,t);
	return sfr;
}
template<class T>
double PezzulliInflowRadialFlow_rSFR<T>::KSCoeff(double sfr){
	return KSN;
}
template<class T>
double PezzulliInflowRadialFlow_rSFR<T>::sigmagas(double R, double t, Grid* rSFR){
	double sfr = find_SFR(R,t,rSFR);
	return pow(sfr/A,1./KSCoeff(sfr));
}
template<class T>
double PezzulliInflowRadialFlow_rSFR<T>::sigmaeffdot(double R, double t, Grid* rSFR){
	double sfr = find_SFR(R,t,rSFR);
	double dsfr = rSFR->t_gradient(R,t);
	double KS = KSCoeff(sfr);
	if(sfr*(1.+pow(sfr,1./KS-2.)*dsfr/pow(A,1./KS)/KS)<0.) return 0.;
	return sfr*(1.+pow(sfr,1./KS-2.)*dsfr/pow(A,1./KS)/KS);
}
template<class T>
double mu_P_integrand(double R, void *p){
	mu_rSFR_st<T> *P=(mu_rSFR_st<T> *) p;
	R = exp(R);
    return pow(R,2.+1./P->P->alpha())*P->P->sigmaeffdot(R,P->t,P->rSFR);
}
template<class T>
double PezzulliInflowRadialFlow_rSFR<T>::mu(double R, double t, Grid* rSFR){
	// radial mass flow rate in units M_sun Gyr^-1
	if(R==0.) return 0.;
	integrator GL(5000);
	// GaussLegendreIntegrator GL2(50);
    mu_rSFR_st<T> Q = {this,t,rSFR};
    // std::cout<<R<<" "<<t<<" "<<" "<<-2.*PI*pow(R,-1./Alpha)*GL.integrate(&mu_P_integrand<T>,log(0.01),log(R),1e-4,&Q)<<" "<<-2.*PI*pow(R,-1./Alpha)*GL2.integrate(&mu_P_integrand<T>,log(0.01),log(R),&Q)<<std::endl;
    return -2.*PI*pow(R,-1./Alpha)*GL.integrate(&mu_P_integrand<T>,log(0.01),log(R),1e-4,&Q);
}
template<class T>
double PezzulliInflowRadialFlow_rSFR<T>::dmudR(double R, double t, Grid* rSFR){
	if(R==0.) return 0.;
	return -mu(R,t,rSFR)/(Alpha*R)-2.*PI*R*sigmaeffdot(R,t,rSFR);
}
template<class T>
double PezzulliInflowRadialFlow_rSFR<T>::acc_rate(double R, double t, Grid* rSFR){
	// surface accretion rate in units M_sun kpc^-2 Gyr^-1
	return -mu(R,t,rSFR)/(2.*PI*Alpha*R*R);
}
template<class T>
double PezzulliInflowRadialFlow_rSFR<T>::flow_rate(double R, double t, Grid* rSFR){
	// radial flow rate in units kpc Gyr^-1
	if(R==0.) return 0.;
	return mu(R,t,rSFR)/sigmagas(R,t,rSFR)/(2.*PI*R);
}
template<class T>
double PezzulliInflowRadialFlow_rSFR<T>::radial_flow_rate(double R, double t, Grid* rSFR){
	// radial flow rate in units M_sun kpc^-2 Gyr^-1
	if(R==0.) return 0.;
	return -dmudR(R,t,rSFR)/(2.*PI*R);
}
template<class T>
double PezzulliInflowRadialFlow_rSFR<T>::operator()(double R, double t, Grid* rSFR){
	return acc_rate(R,t,rSFR);
}
template class PezzulliInflowRadialFlow_rSFR<Inflow>;
template class PezzulliInflowRadialFlow_rSFR<RadialFlow>;
//=============================================================================
// Map for creating shared pointer instances of RadialFlow from string of class name
unique_map< RadialFlow,
            ModelParameters,
            double> radialflow_types ={
    {"Linear",&createInstance<RadialFlow,LinearRadialFlow>},
    {"Pezzulli",&createInstance<RadialFlow,PezzulliInflowRadialFlow_rSFR<RadialFlow>>}
};
// Map for creating shared pointer instances of Inflow from string of class name
unique_map< Inflow,
            ModelParameters,
            double> inflow_types ={
    {"Double",&createInstance<Inflow,DoubleInfallInflow>},
    {"Pezzulli",&createInstance<Inflow,PezzulliInflowRadialFlow_rSFR<Inflow>>}
};
// Map for creating shared pointer instances of Outflow from string of class name
unique_map< Outflow,
            ModelParameters> outflow_types ={
    {"SimpleGalacticFountain",&createInstance<Outflow,SimpleGalacticFountain>}
};
//=============================================================================
