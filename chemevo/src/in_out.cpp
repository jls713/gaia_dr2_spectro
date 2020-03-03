#include "in_out.h"
//=============================================================================
DoubleInfallInflow::DoubleInfallInflow(ModelParameters M,double present_rSFR){
	tf = extract_param(M.parameters["flows"]["inflow"],"FastTimeScale",0.5);
	ts = extract_param(M.parameters["flows"]["inflow"],"SlowTimeScale",8.);
	weight = extract_param(M.parameters["flows"]["inflow"],"Weight",1.);
	Rd = extract_param(M.parameters["fundamentals"],"GasScaleLength",5.);
	pif=1.;
	double pif0 = (*this)(M.parameters["fundamentals"]["SolarRadius"],
	                      M.parameters["fundamentals"]["GalaxyAge"]);
	pif = (double)M.parameters["flows"]["inflow"]["PresentInfallRate"]/pif0;
}
double DoubleInfallInflow::operator()(double R, double t, Grid *rSFR){
	return (exp(-t/ts)*(1.-weight)+exp(-t/tf)*weight)*pif*exp(-R/Rd);
}
//=============================================================================
SimpleGalacticFountain::SimpleGalacticFountain(ModelParameters M){
    auto F = M.parameters["flows"]["outflow"];
    std::vector<std::string> vars = {"Feject_in","Feject_out","TransitionRadius"};
    VecDoub defaults = {0.,0.,3.};
    auto params = extract_params(F, vars, defaults);
    feject_in=params[0];
    feject_out=params[1];
    transR = params[2];
    A=.5*(feject_out+feject_in);
    B=.5*(feject_out-feject_in);
    dR = 0.2;
}
double SimpleGalacticFountain::operator()(double R, double t, double SFR, double GasReturn){
	return (A+B*tanh((R-transR)/dR))*GasReturn;
}
EntrainedWind::EntrainedWind(ModelParameters M){
    eta_wind = extract_param(M.parameters["flows"]["outflow"], "eta_wind", 0.);
}
double EntrainedWind::operator()(double R, double t, double SFR, double GasReturn){
    return eta_wind*SFR;
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
	VGrad = extract_param(M.parameters["flows"]["radialflow"],"Gradient",0.);
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
	Alpha = extract_param(M.parameters["flows"]["radialflow"],
	                      "PezzulliAlpha",0.);
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
GasDumpSimple::GasDumpSimple(ModelParameters M){
	check_param_given(M.parameters["flows"],"gasdump");
	std::vector<std::string> var = {
		"surfacedensity","time","central_radius",
		"radial_width","metallicity","alpha"
	};
	VecDoub defaults = {0., 1., 8.3, 1., -1., 0.};
	auto dd = extract_params(M.parameters["flows"]["gasdump"],var,defaults);
	surfacedensity=dd[0];time=dd[1];central_radius=dd[2];
	radial_width=dd[3];metallicity=dd[4];alpha=dd[5];

	// Put the time near a gridpoint
	double mint=100000., tt;
	auto gas_mass = make_unique<Grid>(M);
	for(int i=0;i<gas_mass->grid_time().size();++i){
		if(fabs(gas_mass->grid_time()[i]-time)<mint){
			tt = gas_mass->grid_time()[i];
			mint = fabs(tt-time);
		}
	}
	time=tt;
}
double GasDumpSimple::operator()(double R, double t, double dt){
	if (abs(t-time)<1e-8)
	    return surfacedensity/dt*exp(-.5*pow((R-central_radius)/radial_width,2.));
	else
		return 0.;
}

double GasDumpSimple::elements(Element E, double R, double t, double dt,
                               SolarAbundances solar){
	if (abs(t-time)<1e-8){
		double metal = solar.Z()*pow(10.,metallicity);
		// Fill each element accounting for initial non-zero alpha enrichment
		auto initial_abundance=solar.scaled_solar_mass_frac(E,metal);
		if(is_alpha_element[E])
			initial_abundance*=pow(10.,alpha);

	    return (*this)(R, t, dt) * initial_abundance;
	}
	else
		return 0.;
}
//=============================================================================
// Map for creating shared pointer instances of RadialFlow from string of class name
unique_map< RadialFlow,
            ModelParameters,
            double> radialflow_types ={
    {"None",&createInstance<RadialFlow,RadialFlowNone>},
    {"Linear",&createInstance<RadialFlow,LinearRadialFlow>},
    {"Pezzulli",&createInstance<RadialFlow,PezzulliInflowRadialFlow_rSFR<RadialFlow>>}
};
// Map for creating shared pointer instances of Inflow from string of class name
unique_map< Inflow,
            ModelParameters,
            double> inflow_types ={
    {"None",&createInstance<Inflow,InflowNone>},
    {"Double",&createInstance<Inflow,DoubleInfallInflow>},
    {"Pezzulli",&createInstance<Inflow,PezzulliInflowRadialFlow_rSFR<Inflow>>}
};
// Map for creating shared pointer instances of Outflow from string of class name
unique_map< Outflow,
            ModelParameters> outflow_types ={
    {"None",&createInstance<Outflow,OutflowNone>},
    {"SimpleGalacticFountain",&createInstance<Outflow,SimpleGalacticFountain>},
    {"EntrainedWind",&createInstance<Outflow,EntrainedWind>}
};
// Map for creating shared pointer instances of GasDump from string of class name
unique_map< GasDump,
            ModelParameters> gasdump_types ={
    {"None",&createInstance<GasDump,GasDumpNone>},
    {"SimpleGasDump",&createInstance<GasDump,GasDumpSimple>}
};
//=============================================================================
