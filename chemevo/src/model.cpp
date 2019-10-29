#include "model.h"
using namespace H5;
unsigned iterMAX_integral=50000;
//=============================================================================
// Setup model
void Model::setup(void){
	//=========================================================================
	LOG(INFO)<<"Setting up model with parameters:"<<params.parameters.dump();
	//=========================================================================
	// 0. Check if parameter values are valid
	if(check_parameters())
		throw std::invalid_argument("Invalid parameters in parameter file\n");
	//=========================================================================
	// 1. Initialise Initial Mass Function, Star Formation Rate,
	//    Stellar Ages, Type Ia Rates and Radial Migration.
	imf = imf_types[params.parameters["fundamentals"]["IMF"]](params);
	sfr = sfr_types[params.parameters["fundamentals"]["SFR"]](params);
	lifetime = life_types[params.parameters["fundamentals"]["lifetimes"]](params);
	rad_mig = rm_types[params.parameters["migration"]["Form"]](params);
	migration = !(params.parameters["migration"]["Form"]=="None");
	typeIaRate = tia_types[params.parameters["typeIa"]["Form"]](params,imf,sfr,rad_mig,lifetime);
	//=========================================================================
	// 2. Initialize yields
	yields = make_unique<YieldsSet>(params);
	agb_yields = !(params.parameters["yields"]["AGB"]=="None");
	typeII_yields = !(params.parameters["yields"]["typeII"]=="None");
	typeIa_yields = !(params.parameters["yields"]["typeIa"]=="None");
	//=========================================================================
	// 3. Initialise Grids of metallicity, H, He and elements
	gas_mass = make_unique<Grid>(params);
	stellar_mass = make_unique<Grid>(params);
	reducedSFR = make_unique<Grid>(params);
	metallicity = make_unique<Grid>(params);

	elements.insert(std::make_pair(0,"H"));
	elements_r.insert(std::make_pair("H",0));
	mass_fraction.push_back(Grid(params));
	elements.insert(std::make_pair(1,"He"));
	elements_r.insert(std::make_pair("He",1));
	mass_fraction.push_back(Grid(params));

	auto i=2;
	for(std::string e:params.parameters["elements"]){
		elements.insert(std::make_pair(i,e));
		elements_r.insert(std::make_pair(e,i));
		++i;
		mass_fraction.push_back(Grid(params));
	}
	//=========================================================================
	// 4. Initialise flows
	outflow = outflow_types[params.parameters["flows"]["outflow"]["Form"]](params);
	double Rs = params.parameters["fundamentals"]["SolarRadius"];
	double ts = params.parameters["fundamentals"]["GalaxyAge"];
	double prSFR=SFR(Rs,ts)-(1.-OutflowFraction(Rs,ts))*GasReturnRate(Rs,ts);
	inflow = inflow_types[params.parameters["flows"]["inflow"]["Form"]](params,prSFR);
	radialflow = radialflow_types[params.parameters["flows"]["radialflow"]["Form"]](params,prSFR);
}
void Model::fill_initial_grids(void){
	double Zinit = params.parameters["fundamentals"]["InitialMetallicity"];
	double AlphaInit = params.parameters["fundamentals"]["InitialAlpha"];
	Zinit *= solar.Z();
        metallicity->set_fixed_t_const(Zinit,0);
	for(auto i=0u;i<mass_fraction.size();++i){
		auto initial_abundance = solar.scaled_solar_mass_frac(elements[i],Zinit);
		if(is_alpha_element[elements[i]])
			initial_abundance*=pow(10.,AlphaInit);
		mass_fraction[i].set_fixed_t_const(initial_abundance,0);
	}
	// Fill gas and star initial grid
	VecDoub gas_radial_dist=gas_mass->grid_radial();
	VecDoub star_radial_dist=gas_mass->grid_radial();

	double Rs = params.parameters["fundamentals"]["SolarRadius"];
	double ts = params.parameters["fundamentals"]["GalaxyAge"];
	double KSN = params.parameters["fundamentals"]["Kennicutt-Schmidt_Coeff"];
	auto F = params.parameters["fundamentals"];
        double A;
	if (F.find("Kennicutt-Schmidt_A") != F.end())
		A = F["Kennicutt-Schmidt_A"];
	else{
		double PresentGasDensity = params.parameters["fundamentals"]["PresentGasDensitySun"];
		double rSFR=SFR(Rs,ts)-(1.-OutflowFraction(Rs,ts))*GasReturnRate(Rs,ts);
		if(rSFR<0.)
			throw std::runtime_error("Reduced SFR<0 at R0 at final time\n");

		A=rSFR/pow(PresentGasDensity,KSN);
	}
	for(auto i=0u;i<gas_radial_dist.size();++i){
		star_radial_dist[i]=SFR(gas_mass->grid_radial()[i],0.);
		gas_radial_dist[i]=pow(star_radial_dist[i]/A,1./KSN);
	}
	reducedSFR->set_fixed_t(star_radial_dist,0);
	gas_mass->set_fixed_t(gas_radial_dist,0);
        LOG(INFO)<<"Grids filled\n";

        if(gasdump){
		double mint=100000., tt;
                for(int i=0;i<gas_mass->grid_time().size();++i)
			if(fabs(gas_mass->grid_time()[i]-gasdumptime)<mint){
				tt = gas_mass->grid_time()[i];
				mint = fabs(gas_mass->grid_time()[i]-gasdumptime);
			}
		gasdumptime=tt;
	}
}

int Model::check_parameters(void){
        if(params.parameters["grids"]["RadialGridPoints"]==1){
            single_zone=true;
            params.parameters["fundamentals"]["SolarRadius"]=1.;
            params.parameters["grids"]["MinimumRadius"]=0.;
            params.parameters["grids"]["MaximumRadius"]=2.;
		    params.parameters["migration"]["Form"]="None";
		    params.parameters["flows"]["radialflow"]["Form"]="None";
		    params.parameters["fundamentals"]["GasScaleLength"]=1;
        }
	auto F = params.parameters["fundamentals"];
	if (F.find("IMF") == F.end()) {
		LOG(INFO)<<"No IMF found in parameters file\n";
		return 1;
	}
	if(imf_types.count(params.parameters["fundamentals"]["IMF"])!=1){
		std::string types = "";
		for(auto i:imf_types) types+=i.first+",";
		LOG(INFO)<<"Invalid IMF:"<<params.parameters["fundamentals"]["IMF"]<<std::endl<<"Options are:"<<types;
		return 1;
	}
	if (F.find("SFR") == F.end()) {
		LOG(INFO)<<"No SFR found in parameters file\n";
		return 1;
	}
	if(sfr_types.count(params.parameters["fundamentals"]["SFR"])!=1){
		std::string types = "";
		for(auto i:sfr_types) types+=i.first+",";
		LOG(INFO)<<"Invalid SFR:"<<params.parameters["fundamentals"]["SFR"]<<std::endl<<"Options are:"<<types;
		return 1;
	}
	if (F.find("lifetimes") == F.end()) {
		LOG(INFO)<<"No lifetimes found in parameters file\n";
		return 1;
	}
	if(life_types.count(params.parameters["fundamentals"]["lifetimes"])!=1){
		std::string types = "";
		for(auto i:life_types) types+=i.first+",";
		LOG(INFO)<<"Invalid Lifetime:"<<params.parameters["fundamentals"]["lifetimes"]<<std::endl<<"Options are:"<<types;
		return 1;
	}
	F = params.parameters;
	if (F.find("typeIa") == F.end()) {
		LOG(INFO)<<"No typeIa found in parameters file\n";
		return 1;
	}
	if(tia_types.count(params.parameters["typeIa"]["Form"])!=1){
		std::string types = "";
		for(auto i:tia_types) types+=i.first+",";
		LOG(INFO)<<"Invalid Type Ia rate:"<<params.parameters["typeIa"]["Form"]<<std::endl<<"Options are:"<<types;
		return 1;
	}
        if(!single_zone){
	if (F.find("migration") == F.end()) {
		LOG(INFO)<<"No migration found in parameters file\n";
		return 1;
	}
	if(rm_types.count(params.parameters["migration"]["Form"])!=1){
		std::string types = "";
		for(auto i:rm_types) types+=i.first+",";
		LOG(INFO)<<"Invalid radial migration type:"<<params.parameters["migration"]["Form"]<<std::endl<<"Options are:"<<types;
		return 1;
	}
	}
	if(inflow_types.count(params.parameters["flows"]["inflow"]["Form"])!=1){
		std::string types = "";
		for(auto i:inflow_types) types+=i.first+",";
		LOG(INFO)<<"Invalid Inflow:"<<params.parameters["flows"]["inflow"]["Form"]<<std::endl<<"Options are:"<<types;
		return 1;
	}
	if(outflow_types.count(params.parameters["flows"]["outflow"]["Form"])!=1){
		std::string types = "";
		for(auto i:outflow_types) types+=i.first+",";
		LOG(INFO)<<"Invalid Outflow:"<<params.parameters["flows"]["outflow"]["Form"]<<std::endl<<"Options are:"<<types;
		return 1;
	}
    if(!single_zone)
	if(radialflow_types.count(params.parameters["flows"]["radialflow"]["Form"])!=1){
		std::string types = "";
		for(auto i:radialflow_types) types+=i.first+",";
		LOG(INFO)<<"Invalid Radial Flow:"<<params.parameters["flows"]["radialflow"]["Form"]<<std::endl<<"Options are:"<<types;
		return 1;
	}
	F = params.parameters["elements"];
	for(auto f: F)
		if(element_index.find(f)==element_index.end()){
			LOG(INFO)<<"Element "<<f<<" not valid."<<std::endl;
		}
	F = params.parameters["flows"];
    if (F.find("gasdump") != F.end()){
		gasdump=true;
		gasdumptime=F["gasdump"]["time"];
		gasdumpsurfacedensity=F["gasdump"]["surfacedensity"];
		gasdumpAlpha=F["gasdump"]["alpha"];
		gasdumpMetal=F["gasdump"]["metallicity"];
	}
	return 0;
}

//=============================================================================
// Run models
int Model::step(unsigned nt, double dt){

	int err = 0;

	// Key Loop
	auto t = gas_mass->grid_time()[nt];
	unsigned nR0=0,nt0=0;
	auto NR = gas_mass->grid_radial().size();
	auto tp = t-.5*dt;

	// First we compute SFR and gas return rate -- this allows us to find the
	// gas needed in inflows in this time step
	#pragma omp parallel for schedule(dynamic)
	for(auto nR=0u;nR<NR;++nR){
		double R = gas_mass->grid_radial()[nR];
		// we set the metallicity to that of the previous time step so small
		// age interpolation works (i.e. doesn't use zero)
		metallicity->set(Z(R,t-dt),nR,nt);
		if(Z(R,t)<0.){
			err=1;
			std::cerr<<"Returning here: "<<" "<<nR<<" ";
			std::cerr<<Z(gas_mass->grid_radial()[nR],t)<<" ";
			std::cerr<<Z(gas_mass->grid_radial()[nR],t-dt)<<std::endl;
			continue;
		}
		double starformrate = SFR(R,t);
		double gas_return = (1.-OutflowFraction(R,t))*GasReturnRate(R,t);
		reducedSFR->set(starformrate-gas_return,nR,nt);
		// if(nR==NR-1) reducedSFR->set(0.,nR,nt);
	}
	if(err){return err;}
        #pragma omp parallel for schedule(dynamic)
	for(auto nR=0u;nR<NR;++nR){
		auto Xi=0.,R=0.,Rdown=0.,Rup=0.;
		auto starformrate=0.,inflowrate=0.,rad_flow_dm=0.,gas_return=0.;
		auto dmdt=0.,dmsdt=0.;
		auto sm_prev=0.,sm=0.;
		auto gm_prev=0.,gm=0.;
		auto rmlower=0.,rmupper=0.,rmlower_1=0.,rmupper_1=0.,gmprev_d=0.;
		auto area=0., area_up=0., area_down=0.;

		R = gas_mass->grid_radial()[nR];

		gm_prev = (*gas_mass)(nR,nt-1);
		starformrate=SFR(R,tp);
		gas_return=starformrate-(*reducedSFR)(R,tp);
	        inflowrate = InflowRate(R,tp);
		//double A = params.parameters["fundamentals"]["Kennicutt-Schmidt_A"];
		//double K = params.parameters["fundamentals"]["Kennicutt-Schmidt_Coeff"];

		//inflowrate = (pow(starformrate/A,K)-gm_prev-(-starformrate+gas_return)*dt)/dt;

		auto outer_gas_mass=0.;
		if(!single_zone){
		Rdown=gas_mass->Rdown(nR);
		Rup=gas_mass->Rup(nR);

		if(nR<NR-1)
			outer_gas_mass=(*gas_mass)(nR+1,nt-1);
		else
			outer_gas_mass=gas_mass->log_extrapolate_high(Rup,nt-1);
			// outer_gas_mass=exp(log((*gas_mass)(NR-1,nt-1))+(log((*gas_mass)(NR-1,nt-1))-log((*gas_mass)(NR-2,nt-1)))/(gas_mass->grid_radial()[NR-1]-gas_mass->grid_radial()[NR-2])*(Rup-gas_mass->grid_radial()[NR-1]));

		rad_flow_dm=RadialFlowRate(gm_prev,outer_gas_mass,R,Rdown,Rup,tp,dt,&err);
		}
        dmdt += -starformrate;
		dmdt += gas_return;
		dmdt += inflowrate;
		dmdt += rad_flow_dm;
		// if(nR==NR-1) dmdt=inflowrate+rad_flow_dm;
		if(((nR>NR-15) or (single_zone&&(nR==NR-1))) and fabs(gas_mass->grid_time()[nt]-gasdumptime)<1e-8 and gasdump){
			double mid,std=1.;
			if(single_zone) mid=gas_mass->grid_radial()[nR];
			else{mid=gas_mass->grid_radial()[NR-7]; std=gas_mass->grid_radial()[NR-5]-mid;}
			dmdt+=gasdumpsurfacedensity/dt*exp(-.5*pow((gas_mass->grid_radial()[nR]-mid)/std,2.));
		}
		gm = gm_prev+dmdt*dt;
                if(gm<0.){starformrate=0.;rad_flow_dm=0.;dmdt=gas_return+inflowrate;gm=gm_prev+dmdt*dt;}

		if(migration)
			gm += rad_mig->convolve(gas_mass.get(),nR,nt)-gm_prev;
		/*if(gm<0.)
			throw std::runtime_error(
			           "More stars being formed than gas available in ring R="
			           +std::to_string(R));
		*/

		// if(migration){
		// 	area_up = gas_mass->annulus_area(nR+1);
		// 	area = gas_mass->annulus_area(nR);
		// 	area_down = gas_mass->annulus_area(nR-1);
		// 	rmlower=RadialMigrationKernel(Rdown,R,dt); // leaving gas moving in
		// 	rmupper=RadialMigrationKernel(Rup,R,dt); // leaving gas moving out
		// 	rmlower_1=RadialMigrationKernel(R,Rdown,dt)*area_down/area; // entering gas moving in
		// 	rmupper_1=RadialMigrationKernel(R,Rup,dt)*area_up/area; // entering gas moving out
		// 	if(nR==0)
		// 		gmprev_d=gas_mass->log_extrapolate_low(Rdown,nt-1);
		// 		// gmprev_d=exp(log((*gas_mass)(0,nt-1))+(log((*gas_mass)(1,nt-1))-log((*gas_mass)(0,nt-1)))/(gas_mass->grid_radial()[1]-gas_mass->grid_radial()[0])*(Rdown-gas_mass->grid_radial()[0]));
		// 	else
		// 		gmprev_d=(*gas_mass)(nR-1,nt-1);
		// 	// if(nR==0)
		// 	// 	gm+=outer_gas_mass*rmupper_1-gm_prev*rmupper;
		// 	// else
		// 		if(nR>0)
		// 		gm+=gmprev_d*rmlower_1+outer_gas_mass*rmupper_1-gm_prev*(rmupper+rmlower);
		// 	// std::cout<<R<<" "<<t<<" "<<rmlower<<" "<<rmupper<<" "<<rmlower_1<<" "<<rmupper_1<<" "<<" "<<Rdown<<" "<<Rup<<" "<<gm_prev<<" "<<gmprev_d<<" "<<gm_prev+dmdt*dt<<" "<<gm<<std::endl;
		// 		std::cout<<rmlower<<" "<<rmupper<<" "<<rmlower_1<<" "<<rmupper_1<<" "<<gmprev_d<<" "<<outer_gas_mass<<" "<<gm_prev<<" "<<Rdown<<" "<<Rup<<" "<<gas_mass->annulus_area(nR)<<" ";
		// }

		// std::cerr<<R<<" "<<t<<" "<<gm<<" "<<starformrate<<" "<<gas_return<<" "<<inflowrate<<" ";
		// std::cerr<<rad_flow_dm<<" "<<reducedSFR->t_gradient(R,tp)<<" "<<outer_gas_mass<<" ";
		//std::cerr<<radialflow->gamma_g(R,Rdown,Rup,tp,dt,reducedSFR.get())<<" "<<radialflow->beta_g(R,Rdown,Rup,tp,dt,reducedSFR.get());
		//std::cerr<<" "<<radialflow->flow_rate((R+Rdown)*.5,tp,reducedSFR.get())<<" "<<radialflow->flow_rate((R+Rup)*.5,tp,reducedSFR.get());

		sm_prev = (*stellar_mass)(nR,nt-1);
		dmsdt = starformrate;
		dmsdt -= gas_return/(1.-OutflowFraction(R,tp));
		sm = sm_prev+dmsdt*dt;

		if(migration and nt>1)
			sm += rad_mig->convolve(stellar_mass.get(),nR,nt)-sm_prev;

		/*if(gm<0.)
			throw std::runtime_error(
			           "More stars being formed than gas available in ring R="
			           +std::to_string(R));
		*/

		gas_mass->set(gm,nR,nt);
		stellar_mass->set(sm,nR,nt);

		for(auto e: elements){
			Xi = mass_fraction[e.first](nR,nt-1);
			dmdt = 0.;
			dmdt = -starformrate*Xi;			  // star formation
			dmdt += (1.-OutflowFraction(R,tp))*EnrichmentRate(e.second,R,tp); // re-enrich
			dmdt += inflowrate*mass_fraction[e.first](nR0,nt0); // inflow
			// for outermost ring we use the initial mass fractions
			if(!single_zone){
			auto e_outer_gas_mass=outer_gas_mass;
			if(nR<NR-1)
				e_outer_gas_mass *= mass_fraction[e.first](nR+1,nt-1);
			else
				e_outer_gas_mass *= mass_fraction[e.first](nR0,nt0);
			rad_flow_dm=RadialFlowRate(gm_prev*Xi,e_outer_gas_mass,R,Rdown,Rup,tp,dt,&err);
			dmdt += rad_flow_dm;
			}

			if(((nR>NR-15) or (single_zone&&(nR==NR-1))) and fabs(gas_mass->grid_time()[nt]-gasdumptime)<1e-8 and gasdump){
        			double metal = solar.Z()*pow(10.,gasdumpMetal);
        			// Fill each element accounting for initial non-zero alpha enrichment
                 		auto initial_abundance=solar.scaled_solar_mass_frac(e.second,metal);
                                if(is_alpha_element[e.second])
                                    initial_abundance*=pow(10.,gasdumpAlpha);
                        	double mid,std=1.;
				if(single_zone) mid=gas_mass->grid_radial()[nR];
				else{mid=gas_mass->grid_radial()[NR-7]; std=gas_mass->grid_radial()[NR-5]-mid;}
				dmdt+=gasdumpsurfacedensity*initial_abundance/dt*exp(-.5*pow((gas_mass->grid_radial()[nR]-mid)/std,2.));
        		}
			// else dmdt=inflowrate*mass_fraction[e.first](nR0,nt0);
			// updated mass of element e
			// if(nR==NR-1) mass_fraction[e.first].set(Xi,nR,nt);

			// double mig=0., lowerXi=0.;
			// if(migration){
			// 	if(nR==0)
			// 		lowerXi=mass_fraction[e.first](0,nt-1);
			// 	else
			// 		lowerXi=mass_fraction[e.first](nR-1,nt-1);
			// 	if(nR>0)
			// 	mig=gmprev_d*lowerXi*rmlower_1+e_outer_gas_mass*rmupper_1-gm_prev*(rmupper+rmlower)*Xi;
			// 	// if(nR==0)
			// 	// 	mig=e_outer_gas_mass*rmupper_1-gm_prev*rmupper*Xi;
			// }

			double mass_now = Xi*gm_prev+dmdt*dt;
			if(migration)
				mass_now+=rad_mig->convolve_massfrac(gas_mass.get(),&mass_fraction[e.first],nR,nt)-Xi*gm_prev;

			mass_fraction[e.first].set(mass_now/gm,nR,nt);
		}
		metallicity->set(1.-X(R,t)-Y(R,t),nR,nt);
		if(Z(R,t)<0.) {
			metallicity->set(Z(R,t-dt),nR,nt);
			mass_fraction[element_index["H"]].set(
			              mass_fraction[element_index["H"]](nR,nt-1),nR,nt);
			mass_fraction[element_index["He"]].set(
			              mass_fraction[element_index["He"]](nR,nt-1),nR,nt);
			err=1;
			std::cerr<<"Problem "<<Z(R,t-dt)<<" "<<Z(R,t)<<std::endl;
		}
		// std::cout<<R<<" "<<t<<" "<<gm<<" "<<Z(R,t)<<std::endl;
		// std::cerr<<" "<<X(R,t)<<" "<<Y(R,t)<<" "<<Z(R,t)<<std::endl;
	}
	if(err) return err;
	int iteratemax=1; auto gm_it=gas_mass->grid_fixed_t(nt);
	std::vector<int> not_done(gm_it.size(),true);
	for(int N=0;N<iteratemax;++N){
		// We recompute the reducedSFR as the metallicity at the current time
		// is slightly different
		#pragma omp parallel for schedule(dynamic)
		for(auto nR=0u;nR<NR;++nR){
			if(not_done[nR]==false) continue;
			double R = gas_mass->grid_radial()[nR];
			double starformrate = SFR(R,t);
			double gas_return = (1.-OutflowFraction(R,t))*GasReturnRate(R,t);
			reducedSFR->set(starformrate-gas_return,nR,nt);
		}
		#pragma omp parallel for schedule(dynamic)
		for(auto nR=0u;nR<NR;++nR){
			if(not_done[nR]==false) continue;
			auto Xi=0.,R=0.,Rdown=0.,Rup=0.,Xihere=0.;
			auto starformrate=0.,inflowrate=0.,rad_flow_dm=0.,gas_return=0.;
			auto dmdt=0.,dmsdt=0.;
			auto sm_prev=0.,sm=0.;
			auto gm_prev=0.,gm=0.,gmhere=0.;
			auto rmlower=0.,rmupper=0.,rmlower_1=0.,rmupper_1=0.,gmprev_d=0.;
			auto area=0., area_up=0., area_down=0.;

			R = gas_mass->grid_radial()[nR];

			gm_prev = (*gas_mass)(nR,nt-1);
			starformrate=SFR(R,tp);
			gas_return=starformrate-(*reducedSFR)(R,tp);
			inflowrate = InflowRate(R,tp);
		//double A = params.parameters["fundamentals"]["Kennicutt-Schmidt_A"];
		//double K = params.parameters["fundamentals"]["Kennicutt-Schmidt_Coeff"];

		//inflowrate = (pow(starformrate/A,K)-gm_prev-(-starformrate+gas_return)*dt)/dt;
			auto outer_gas_mass=0.;
			if(!single_zone){
			Rdown=gas_mass->Rdown(nR);
			Rup=gas_mass->Rup(nR);

			if(nR<NR-1)
				outer_gas_mass=(*gas_mass)(Rup,tp);
			else
				outer_gas_mass=gas_mass->log_extrapolate_high(Rup,tp);
			gmhere=(*gas_mass)(R,tp);

			rad_flow_dm=RadialFlowRate(gmhere,outer_gas_mass,R,Rdown,Rup,tp,dt,&err);
			}
                        dmdt += -starformrate;
			dmdt += gas_return;
			dmdt += inflowrate;
			dmdt += rad_flow_dm;

			if(((nR>NR-15) or (single_zone&&(nR==NR-1))) and fabs(gas_mass->grid_time()[nt]-gasdumptime)<1e-8 and gasdump){
                        	double mid,std=1.;
				if(single_zone) mid=gas_mass->grid_radial()[nR];
				else{mid=gas_mass->grid_radial()[NR-7]; std=gas_mass->grid_radial()[NR-5]-mid;}
				dmdt+=gasdumpsurfacedensity/dt*exp(-.5*pow((gas_mass->grid_radial()[nR]-mid)/std,2.));
			}
			gm = gm_prev+dmdt*dt;

                        if(gm<0.){starformrate=0.;rad_flow_dm=0.;dmdt=gas_return+inflowrate;gm=gm_prev+dmdt*dt;}

			// if(migration){
			// 	area_up = gas_mass->annulus_area(nR+1);
			// 	area = gas_mass->annulus_area(nR);
			// 	area_down = gas_mass->annulus_area(nR-1);
			// 	rmlower=RadialMigrationKernel(Rdown,R,dt); // leaving gas moving in
			// 	rmupper=RadialMigrationKernel(Rup,R,dt); // leaving gas moving out
			// 	rmlower_1=RadialMigrationKernel(R,Rdown,dt)*area_down/area; // entering gas moving in
			// 	rmupper_1=RadialMigrationKernel(R,Rup,dt)*area_up/area; // entering gas moving out
			// 	if(nR==0)
			// 		gmprev_d=gas_mass->log_extrapolate_low(Rdown,tp);
			// 	else
			// 		gmprev_d=(*gas_mass)(Rdown,tp);
			// 	// if(nR==0)
			// 	// 	gm+=outer_gas_mass*rmupper_1-gm_prev*rmupper;
			// 	// else
			// 	if(nR>0)
			// 		gm+=gmprev_d*rmlower_1+outer_gas_mass*rmupper_1-gmhere*(rmupper+rmlower);
			// 	std::cout<<rmlower<<" "<<rmupper<<" "<<rmlower_1<<" "<<rmupper_1<<" "<<gmprev_d<<" "<<outer_gas_mass<<" "<<gm_prev<<" "<<Rdown<<" "<<Rup<<" "<<gas_mass->annulus_area(nR)<<" ";
			// }

			// std::cerr<<R<<" "<<t<<" "<<gm<<" "<<starformrate<<" "<<gas_return<<" "<<inflowrate<<" ";
			// std::cerr<<rad_flow_dm<<" "<<reducedSFR->t_gradient(R,tp)<<" "<<outer_gas_mass<<" ";
			//std::cerr<<radialflow->gamma_g(R,Rdown,Rup,tp,dt,reducedSFR.get())<<" "<<radialflow->beta_g(R,Rdown,Rup,tp,dt,reducedSFR.get())<<" ";
			//std::cerr<<radialflow->flow_rate((R+Rdown)*.5,tp,reducedSFR.get())<<" "<<radialflow->flow_rate((R+Rup)*.5,tp,reducedSFR.get());

			sm_prev = (*stellar_mass)(nR,nt-1);
			dmsdt = starformrate;
			dmsdt -= gas_return/(1.-OutflowFraction(R,tp));
			sm = sm_prev+dmsdt*dt;

			if(migration){
				gm += rad_mig->convolve(gas_mass.get(),nR,nt)-gm_prev;
				// std::cout<<gm-dmdt*dt-gm_prev<<" ";
				if(nt>1) sm += rad_mig->convolve(stellar_mass.get(),nR,nt)-sm_prev;
			}

			/*if(gm<0.)
				throw std::runtime_error(
				           "More stars being formed than gas available in ring R="
				           +std::to_string(R));
			*/
		        //if(gm<0.) gm=0.;
			gas_mass->set(gm,nR,nt);
			stellar_mass->set(sm,nR,nt);

			for(auto e: elements){
				Xi = mass_fraction[e.first](nR,nt-1);
				Xihere = mass_fraction[e.first](R,tp);
				dmdt = 0.;
				dmdt = -starformrate*Xi;			  // star formation
				dmdt += (1.-OutflowFraction(R,tp))*EnrichmentRate(e.second,R,tp); // re-enrich
				dmdt += inflowrate*mass_fraction[e.first](nR0,nt0); // inflow
				// for outermost ring we use the initial mass fractions
                         	if(!single_zone){
			        auto e_outer_gas_mass=outer_gas_mass;
				if(nR<NR-1)
					e_outer_gas_mass *= mass_fraction[e.first](Rup,tp);
				else
					e_outer_gas_mass *= mass_fraction[e.first](nR0,nt0);
				rad_flow_dm=RadialFlowRate(gmhere*Xihere,e_outer_gas_mass,R,Rdown,Rup,tp,dt,&err);
				dmdt += rad_flow_dm;
				}

				if(((nR>NR-15) or (single_zone&&(nR==NR-1))) and fabs(gas_mass->grid_time()[nt]-gasdumptime)<1e-8 and gasdump){
                                	double metal = solar.Z()*pow(10.,gasdumpMetal);
                                	// Fill each element accounting for initial non-zero alpha enrichment
                                	auto initial_abundance=solar.scaled_solar_mass_frac(e.second,metal);
                                	if(is_alpha_element[e.second])
                                		initial_abundance*=pow(10.,gasdumpAlpha);
					double mid,std=1.;
					if(single_zone) mid=gas_mass->grid_radial()[nR];
					else{mid=gas_mass->grid_radial()[NR-7]; std=gas_mass->grid_radial()[NR-5]-mid;}
					dmdt+=gasdumpsurfacedensity*initial_abundance/dt*exp(-.5*pow((gas_mass->grid_radial()[nR]-mid)/std,2.));
				}
				double mass_now = Xi*gm_prev+dmdt*dt;
				if(migration){
					mass_now+=rad_mig->convolve_massfrac(gas_mass.get(),&mass_fraction[e.first],nR,nt)-Xi*gm_prev;
					// std::cout<<mass_now-dmdt*dt-Xi*gm_prev<<" ";
				}
				mass_fraction[e.first].set(mass_now/gm,nR,nt);
			}

			metallicity->set(1.-X(R,t)-Y(R,t),nR,nt);
			if(Z(R,t)<0.) {
				metallicity->set(Z(R,tp),nR,nt); err=1;
			        mass_fraction[element_index["H"]].set(mass_fraction[element_index["H"]](nR,nt-1),nR,nt);
			        mass_fraction[element_index["He"]].set(mass_fraction[element_index["He"]](nR,nt-1),nR,nt);
			}

			// std::cerr<<" "<<X(R,t)<<" "<<Y(R,t)<<" "<<Z(R,t)<<std::endl;
			if(fabs((gm-gm_it[nR])/gm)<0.005) not_done[nR]=0;
			else not_done[nR]=1;
		}
		int done_sum=0;
		for_each(begin(not_done),end(not_done),[&done_sum](int p){done_sum+=p;});
		if(done_sum==0) break;
		else gm_it=gas_mass->grid_fixed_t(nt);
	}
	return err;
}

//=============================================================================
// Run models
void Model::simple_step(unsigned nt, double dt){
	// Key Loop
	auto t = gas_mass->grid_time()[nt];
	unsigned nR0=0,nt0=0;
	auto NR = gas_mass->grid_radial().size();
	auto tp = t-.5*dt;

	for(auto nR=0u;nR<NR;++nR){
		double c = rad_mig->convolve(gas_mass.get(),nR,nt);
		gas_mass->set(c,nR,nt);
		std::cout<<gas_mass->grid_radial()[nR]<<" "<<t<<" "<<(*gas_mass)(nR,nt)<<std::endl;
	}
	// ForwardSolver CN(rad_mig,NR);
	// CN.step(gas_mass.get(),nt,dt);
	// for(auto nR=0u;nR<NR;++nR)
	// 	std::cout<<gas_mass->grid_radial()[nR]<<" "<<t<<" "<<(*gas_mass)(nR,nt)<<std::endl;
/*
	CrankNicolsonSolver CN(rad_mig,NR);
	CN.step(gas_mass.get(),nt,dt);
	for(auto nR=0u;nR<NR;++nR)
		std::cout<<gas_mass->grid_radial()[nR]<<" "<<t<<" "<<(*gas_mass)(nR,nt)<<std::endl;
*/
/*	#pragma omp parallel for schedule(dynamic)
	for(auto nR=0u;nR<NR;++nR){
		auto Xi=0.,R=0.,Rdown=0.,Rup=0.;
		auto starformrate=0.,inflowrate=0.,rad_flow_dm=0.,gas_return=0.;
		auto dmdt=0.,dmsdt=0.;
		auto sm_prev=0.,sm=0.;
		auto gm_prev=0.,gm=0.;
		auto rmlower=0.,rmupper=0.,rmlower_1=0.,rmupper_1=0.,gmprev_d=0.;
		double width_up, width_down, width;
		R = gas_mass->grid_radial()[nR];

		gm = (*gas_mass)(nR,nt-1);
		gm_prev = (*gas_mass)(nR,nt-1);

		Rdown=gas_mass->Rdown(nR);
		Rup=gas_mass->Rup(nR);

		auto outer_gas_mass=0.;
		if(nR<NR-1)
			outer_gas_mass=(*gas_mass)(nR+1,nt-1);
		else
			outer_gas_mass=gas_mass->log_extrapolate_high(Rup,nt-1);

		if(migration){
			width_up = gas_mass->annulus_width(nR+1);
			width = gas_mass->annulus_width(nR);
			width_down = gas_mass->annulus_width(nR-1);
			rmlower=RadialMigrationKernel(Rdown,R,dt)*width; // leaving gas moving in
			rmupper=RadialMigrationKernel(Rup,R,dt)*width; // leaving gas moving out
			rmlower_1=RadialMigrationKernel(R,Rdown,dt)*width_down; // entering gas moving in
			rmupper_1=RadialMigrationKernel(R,Rup,dt)*width_up; // entering gas moving out

			if(nR>0){
				gmprev_d=(*gas_mass)(nR-1,nt-1);
				gm+=gmprev_d*rmlower_1+outer_gas_mass*rmupper_1-gm_prev*(rmupper+rmlower);
			}
			else{
				// zero flux boundary condition
				gm+=outer_gas_mass*rmupper_1-gm_prev*rmupper;
			}
			std::cout<<R<<" "<<t<<" "<<gm<<" "<<rmlower<<std::endl;
		}
		gas_mass->set(gm,nR,nt);
	}*/
}

void Model::expand_grids(unsigned nt, double t){
	gas_mass->add_time(nt,t);
	stellar_mass->add_time(nt,t);
	reducedSFR->add_time(nt,t);
	metallicity->add_time(nt,t);
	for(unsigned mfn=0;mfn<mass_fraction.size();++mfn)
		mass_fraction[mfn].add_time(nt,t);
}

void Model::run(void){
	bool logspace=false;
	auto F = params.parameters["grids"];
	if (F.find("LogAgeGrid") != F.end())
		logspace=F["LogAgeGrid"];
	VecDoub times = gas_mass->grid_time();
	auto NR = gas_mass->grid_radial().size();
	for(unsigned t=1;t<times.size();++t){
		int err;
		if(times[t]-times[t-1]>0.2)
			err=1;
		else
			err=step(t,times[t]-times[t-1]);
		// if there is an error, we bisect time step
		if(err==1)
			while(err==1){
				expand_grids(t,.5*(times[t-1]+times[t]));
				times.insert(times.begin()+t,.5*(times[t-1]+times[t]));
				err=step(t,times[t]-times[t-1]);
			}

		// simple_step(t,times[t]-times[t-1]);
		if(logspace)
			printProgBar((int)(100.*(log(times[t])-log(times[1]))/(log(times.back())-log(times[1]))));
		else
			printProgBar((int)(100.*(times[t]/times.back())));

	}
	std::cout<<std::endl;
	// now scale wrt solar
	metallicity->log10scale(solar.Z());
	for(auto e:elements)
		mass_fraction[e.first].log10scale(solar.mass_fraction(e.second));
}

//=============================================================================
// Output model properties in HDF5 format
void Model::write_properties(H5File &fout){
	VecDoub sfr_grids,inflowss,snII_r,snIa_r;
	double Rs = params.parameters["fundamentals"]["SolarRadius"];
        for(auto t:gas_mass->grid_time()){
		sfr_grids.push_back(SFR(Rs,t));
		inflowss.push_back(InflowRate(Rs,t));
		snII_r.push_back(SNIIRate(Rs,t));
		snIa_r.push_back(SNIaRate(Rs,t));
	}
	hdf5_write_1D_vector(fout,sfr_grids,"SFR");
	hdf5_write_1D_vector(fout,inflowss,"Inflow");
	hdf5_write_1D_vector(fout,snII_r,"SNII");
	hdf5_write_1D_vector(fout,snIa_r,"SNIa");
}
void Model::write(std::string filename){
	H5File fout(filename,H5F_ACC_TRUNC);
	params.write_hdf5(fout);
	gas_mass->write_radial_grid_hdf5(fout);
	gas_mass->write_time_grid_hdf5(fout);
	gas_mass->write_hdf5(fout,"Mgas");
	stellar_mass->write_hdf5(fout,"Mstar");
	write_properties(fout);
	metallicity->write_hdf5(fout,"Z");
	for(auto i=0u;i<elements.size();++i)
		mass_fraction[i].write_hdf5(fout,elements[i]);
	return;
}
//=============================================================================
double Model::X(double R, double t){ return mass_fraction[0](R,t);}
double Model::Y(double R, double t){ return mass_fraction[1](R,t);}
double Model::Z(double R, double t){ return (*metallicity)(R,t);}
//=============================================================================
// Simple interface routines to base classes
double Model::SNIaRate(double R, double t){ return (*typeIaRate)(R,t);}
double Model::SFR(double R, double t){ return (*sfr)(R,t);}
double Model::IMF(double m){ return (*imf)(m);}
double Model::InflowRate(double R, double t){return (*inflow)(R,t,reducedSFR.get());}
double Model::OutflowFraction(double R, double t){return (*outflow)(R,t);}
double Model::RadialFlowRate(double m, double mup, double R, double Rdown, double Rup, double t, double dt, int*err){
	// if(Rup>reducedSFR->grid_radial().back()){
	// 	double dR = (R-Rdown);
	// 	double R4 = R-dR*3.;
	// 	double R3 = R-dR*2.;
	// 	double R2 = R-dR;
	// 	// double gamma_g1=radialflow->gamma_g(R3,R4,R2,t,reducedSFR.get());
	// 	// double beta_g1=radialflow->beta_g(R3,R4,R2,t,reducedSFR.get());
	// 	// double gamma_g2=radialflow->gamma_g(R2,R3,R,t,reducedSFR.get());
	// 	// double beta_g2=radialflow->beta_g(R2,R3,R,t,reducedSFR.get());
	// 	// double gb = (gamma_g2-beta_g2)+(gamma_g2-gamma_g1-beta_g2+beta_g1)/(R2-R3)*(R-R2);
	// 	double r1 = radialflow->dMdt((*gas_mass)(R3,t),(*gas_mass)(R2,t),R3,R4,R2,t,reducedSFR.get());
	// 	double r2 = radialflow->dMdt((*gas_mass)(R2,t),(*gas_mass)(R,t),R2,R3,R,t,reducedSFR.get());
	// 	// double g = (gamma_g2)+(gamma_g2-gamma_g1)/(R2-R3)*(R-R2);
	// 	// return gb*m+g*(mup-m);
	// 	return r2+(r2-r1)/(R2-R3)*(R-R2);
	// }
	return radialflow->dMdt(m,mup,R,Rdown,Rup,t,dt,reducedSFR.get(),err);
}
double Model::RadialMigrationKernel(double R, double Rp,double t){return (*rad_mig)(R,Rp,t);}
double Model::EjectedMass(Element E, double M, double Z){
	return yields->mass(E,M,Z);}
double Model::TotalEjectedMass(double M, double Z){
	return yields->mass_ejected(M,Z);}
double Model::Lifetime(double m, double Z){ return (*lifetime)(m,Z);}
double Model::Mass_Star_Dying_Now(double age, double Z){
	return lifetime->mass_star_dying_now(age,Z);
}
double Model::dMdt(double M, double age, double Z){
	return M*lifetime->dlogMdlogt(age,Z)/age;
}

//=============================================================================
double _snIIrate(double m,void *p){
	Rate_st *P = (Rate_st *)p;
	return P->M->IMF(m)*P->M->SFR(P->R,P->t-P->M->Lifetime(m,0.02));
}

double Model::SNIIRate(double R, double t){
	GaussLegendreIntegrator GL(100);
	auto mmin = params.min_typeII_SN_mass();
	auto mmax = params.max_mass();
	auto mm = lifetime->mass_star_dying_now(t,Z(R,t));
	if(mm>mmax) return 0.;
	if(mm>mmin) mmin=mm;
	Rate_st P(this,R,t,"",0.);
	return GL.integrate(&_snIIrate,mmin,mmax,&P);
}
//=============================================================================
// Death Rates
double _death_rate(double tp,void *p){
	Rate_st *P = (Rate_st *)p;
	if(tp>P->t) return 0.;
	auto Z = P->M->Z(P->R,tp);
	auto age = P->t-tp;
	auto M = P->M->Mass_Star_Dying_Now(age,Z);
    if(P->which=="AGB" and M>P->mthresh)
    	return 0.;
   	else if(P->which=="TypeII" and M<P->mthresh)
   		return 0.;
	auto dmdt = P->M->dMdt(M,age,Z);
	return P->M->IMF(M)*P->M->SFR(P->R,tp)*(-dmdt);
}

int _death_rate_2D(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
    double y2[2];
    Rate_st_2D *P = (Rate_st_2D *) fdata;
    for(int i=0;i<2;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
    auto tp = y2[0], Rp = y2[1];
    if(tp>P->t) return 0.;
    auto Z = P->M->Z(P->R,tp);
    auto age = P->t-tp;
    auto M = P->M->Mass_Star_Dying_Now(age,Z);
    if(P->which=="AGB" and M>P->mthresh)
    	return 0.;
   	else if(P->which=="TypeII" and M<P->mthresh)
   		return 0.;
    auto dmdt = P->M->dMdt(M,age,Z);
    fval[0]=P->M->IMF(M)*P->M->SFR(Rp,tp)*(-dmdt)*P->M->RadialMigrationKernel(P->R,Rp,age);
    return 0;
}

double Model::SwitchDeathRate(double R, double t, double tmin, double tmax,std::string which){
	double mthresh = params.parameters["typeII"]["Min_typeII_SN_mass"];
	if(params.parameters["migration"]["Form"]=="None"){
		Rate_st P({this,R,t,which,mthresh});
		//GaussLegendreIntegrator GL(400);
		integrator GL(iterMAX_integral);
		return GL.integrate(&_death_rate,tmin,tmax,1e-4,&P);
	}
	else{
		double Rmin=0.,Rmax=30.;
		Rate_st_2D P(this,R,t,{tmin,Rmin},{tmax,Rmax},which,mthresh);
		double err;
		return integrate(&_death_rate_2D,&P,1e-3,0,"Divonne",&err,
		                 "DeathRate "+which);
	}
}

double Model::TestSNIIRate(double R, double t){
	auto tmin = t-Lifetime(
	                params.parameters["typeII"]["Min_typeII_SN_mass"],Z(R,0.));
	if(tmin<0.) tmin=0.;
	auto tmax = t-Lifetime(params.max_mass(),Z(R,t));
	if(tmax<0.) return 0.;
	return SwitchDeathRate(R,t,tmin,tmax,"SNII");
}
double Model::AGBDeathRate(double R, double t){
	auto tmin = t-Lifetime(params.min_mass(),Z(R,0.));
	if(tmin<0.) tmin=0.;
	auto tmax = t-Lifetime(params.parameters["typeII"]["Min_typeII_SN_mass"],
	                       Z(R,t));
	if(tmax<0.) return 0.;
	return SwitchDeathRate(R,t,tmin,tmax,"AGB");
}
double Model::DeathRate(double R, double t){
	return (agb_yields?AGBDeathRate(R,t):0.)+(typeII_yields?SNIIRate(R,t):0.);
}

//=============================================================================
// Enrichment rates
double _abundance_enrichment_rate(double tp,void *p){
	EnrichRate_st *P = (EnrichRate_st *)p;
	if(tp>P->t) return 0.;
	auto Z = P->M->Z(P->R,tp);
	auto age = P->t-tp;
	auto M = P->M->Mass_Star_Dying_Now(age,Z);
    if(P->which=="AGB" and M>P->mthresh)
    	return 0.;
   	else if(P->which=="TypeII" and M<P->mthresh)
   		return 0.;
	auto dmdt = P->M->dMdt(M,age,Z);
	return P->M->EjectedMass(P->E,M,Z)*P->M->IMF(M)*P->M->SFR(P->R,tp)*(-dmdt);
}
int _abundance_rate_2D(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
    double y2[2];
    EnrichRate_st_2D *P = (EnrichRate_st_2D *) fdata;
    for(int i=0;i<2;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
    double tp = y2[0], Rp = y2[1];
    if(tp>P->t) return 0.;
    double Z = P->M->Z(P->R,tp);
    double age = P->t-tp;
    double M = P->M->Mass_Star_Dying_Now(age,Z);
    if(P->which=="AGB" and M>P->mthresh)
    	return 0.;
   	else if(P->which=="TypeII" and M<P->mthresh)
   		return 0.;
    auto dmdt = P->M->dMdt(M,age,Z);
    fval[0]=P->M->EjectedMass(P->E,M,Z)*P->M->IMF(M)*P->M->SFR(Rp,tp)*(-dmdt)*P->M->RadialMigrationKernel(P->R,Rp,age);
    // if(P->which=="AGB")
    // 	std::cout<<Rp<<" "<<tp<<" "<<P->E<<" "<<fval[0]<<std::endl;
        return 0;
}


double Model::SwitchEnrichmentRate(Element E, double R, double t, double tmin, double tmax,std::string which){
	double mthresh = params.parameters["typeII"]["Min_typeII_SN_mass"];
	// We assume migration is sufficiently weak that TypeIIs do not migrate
	// at all.
	if(params.parameters["migration"]["Form"]=="None" or which=="TypeII"){
		EnrichRate_st P(this,R,t,E,which,mthresh);
		// GaussLegendreIntegrator GL(200);
		// return GL.integrate(&_abundance_enrichment_rate,tmin,tmax,&P);
		integrator GL(iterMAX_integral);
		return GL.integrate(&_abundance_enrichment_rate,tmin,tmax,1e-4,&P);
	}
	else{
		double Rmin=0.,Rmax=30.;
		VecDoub xmin = {tmin,Rmin};
		VecDoub xmax = {tmax,Rmax};
		EnrichRate_st_2D P(this,R,t,E,xmin,xmax,which,mthresh);
		double err;
		return integrate(&_abundance_rate_2D,&P,1e-3,0,"Divonne",&err,
		                 "EnrichmentRate "+which);
	}
}
double Model::TypeIIEnrichmentRate(Element E, double R, double t){
	auto tmin = t-Lifetime(params.parameters["typeII"]["Min_typeII_SN_mass"],
	                       Z(R,0.));
	if(tmin<0.) tmin=0.;
	auto tmax = t-Lifetime(params.max_mass(),Z(R,t));
	if(tmax<0.) return 0.;
	return SwitchEnrichmentRate(E,R,t,tmin,tmax,"TypeII");
}

double Model::AGBEnrichmentRate(Element E, double R, double t){
	auto tmin = t-Lifetime(params.min_mass(),Z(R,0.));
	if(tmin<0.) tmin=0.;
	auto tmax = t-Lifetime(params.parameters["typeII"]["Min_typeII_SN_mass"],
	                       Z(R,t));
	if(tmax<0.) return 0.;
	return SwitchEnrichmentRate(E,R,t,tmin,tmax,"AGB");
}
double Model::TypeIaEnrichmentRate(Element E, double R, double t){
	return yields->typeIa_ejectedmass(E)*SNIaRate(R,t);
}
double Model::EnrichmentRate(Element E, double R, double t){
	return (typeIa_yields?TypeIaEnrichmentRate(E,R,t):0.)
		  +(typeII_yields?TypeIIEnrichmentRate(E,R,t):0.)
		  +(agb_yields?AGBEnrichmentRate(E,R,t):0.);
}

//=============================================================================
// Mass of gas return rate
double _gasreturn_rate(double tp,void *p){
	Rate_st *P = (Rate_st *)p;
	if(tp>P->t) return 0.;
	auto Z = P->M->Z(P->R,tp);
	auto age = P->t-tp;
	auto M = P->M->Mass_Star_Dying_Now(age,Z);
    if(P->which=="AGB" and M>P->mthresh)
    	return 0.;
   	else if(P->which=="TypeII" and M<P->mthresh)
   		return 0.;
	auto dmdt = P->M->dMdt(M,age,Z);
	return P->M->TotalEjectedMass(M,Z)*P->M->IMF(M)*P->M->SFR(P->R,tp)*(-dmdt);
}
int _gasreturn_rate_2D(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
    double y2[2];
    Rate_st_2D *P = (Rate_st_2D *) fdata;
    for(int i=0;i<2;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
    auto tp = y2[0], Rp = y2[1];
    if(tp>P->t) return 0.;
    auto Z = P->M->Z(P->R,tp);
    auto age = P->t-tp;
    auto M = P->M->Mass_Star_Dying_Now(age,Z);
    if(P->which=="AGB" and M>P->mthresh)
    	return 0.;
   	else if(P->which=="TypeII" and M<P->mthresh)
   		return 0.;
    auto dmdt = P->M->dMdt(M,age,Z);
    fval[0]=P->M->TotalEjectedMass(M,Z)*P->M->IMF(M)*P->M->SFR(P->R,tp)*(-dmdt)*P->M->RadialMigrationKernel(P->R,Rp,age);
    return 0;
}

double Model::SwitchGasReturnRate(double R, double t, double tmin, double tmax,std::string which){
	double mthresh = params.parameters["typeII"]["Min_typeII_SN_mass"];
	if(params.parameters["migration"]["Form"]=="None"){
		Rate_st P({this,R,t,which,mthresh});
		//GaussLegendreIntegrator GL(300);
		integrator GL(iterMAX_integral);
		return GL.integrate(&_gasreturn_rate,tmin,tmax,1e-4,&P);
	}
	else{
		double Rmin=0.,Rmax=30.;
		// if(params.parameters["migration"]["Form"]=="Gaussian"){
		// 	Rmin=R-4.*(double)params.parameters["migration"]["sigmaR"];
		// 	Rmin=(Rmin<0.?0.:Rmin);
		// 	Rmax=R+4.*(double)params.parameters["migration"]["sigmaR"];
		// }
		Rate_st_2D P(this,R,t,{tmin,Rmin},{tmax,Rmax},which,mthresh);
		double err;
		return integrate(&_gasreturn_rate_2D,&P,1e-3,0,"Divonne",&err,
		                 "GasReturn "+which);
	}
}
double Model::TypeIIGasReturnRate(double R, double t){
	auto tmin = t-Lifetime(params.parameters["typeII"]["Min_typeII_SN_mass"],
	                       Z(R,0.));
	if(tmin<0.) tmin=0.;
	auto tmax = t-Lifetime(params.max_mass(),Z(R,t));
	if(tmax<0.) return 0.;
	return SwitchGasReturnRate(R,t,tmin,tmax,"Type II");
}

double Model::AGBGasReturnRate(double R, double t){
	auto tmin = t-Lifetime(params.min_mass(),Z(R,0.));
	if(tmin<0.) tmin=0.;
	auto tmax = t-Lifetime(params.parameters["typeII"]["Min_typeII_SN_mass"],
	                       Z(R,t));
	if(tmax<0.) return 0.;
	return SwitchGasReturnRate(R,t,tmin,tmax,"AGB");
}
double Model::TypeIaGasReturnRate(double R, double t){
	return Mch*SNIaRate(R,t);
}
double Model::GasReturnRate(double R, double t){
	return (typeIa_yields?TypeIaGasReturnRate(R,t):0.)
		  +(typeII_yields?TypeIIGasReturnRate(R,t):0.)
		  +(agb_yields?AGBGasReturnRate(R,t):0.);
}

//=============================================================================
//=============================================================================
