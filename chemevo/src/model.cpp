#include "model.h"
using namespace H5;
unsigned iterMAX_integral=50000;
//=============================================================================
// Setup model
void Model::setup(void){
	//=========================================================================
	LOG(INFO)<<"Setting up model with parameters:"<<params.parameters.dump();
	params.pretty_print(std::cout);
	//=========================================================================
	// 0. Check if parameter values are valid
	if(check_parameters())
		throw std::invalid_argument("Invalid parameters in parameter file\n");
	// Initialize iteration number per step
	iteratemax = extract_param(params.parameters["fundamentals"],
	                           "max_iterations",2);
	tol = extract_param(params.parameters["fundamentals"],
	                           "tolerance",0.005);
	std::string default_string = "Asplund";
	std::string solar_string = extract_param(params.parameters["fundamentals"],"solar",default_string);
	params.parameters["fundamentals"]["solar"]=solar_string;
	solar = solar_types[solar_string](params);
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
	logspace=extract_param(params.parameters["grids"],"LogAgeGrid",true);
	params.parameters["grids"]["LogAgeGrid"]=logspace;
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

	// Set up warm phase
	// warm_cold_ratio = f_g in the overleaf
	warm_cold_ratio = extract_param(params.parameters["fundamentals"],
	                                "warm_cold_ratio",0.);
        // warm_cooling_time
	warm_cooling_time = extract_param(params.parameters["fundamentals"],
	                                "warm_cooling_time",0.);
	// outflow_warm_fraction = f_o in the overleaf
	outflow_warm_fraction = extract_param(params.parameters["flows"]["outflow"],
	                                      "warm_fraction_retained",0.);
	if(warm_cold_ratio>1.)
		throw std::invalid_argument("Invalid warm/cold ratio (>1)\n");
	if(warm_cold_ratio>0. or outflow_warm_fraction>0.){
		use_warm_phase=true;
		warm_gas_mass = make_unique<Grid>(params);
		for(auto mf: mass_fraction)
			mass_fraction_warm.push_back(Grid(params));
	}
	else use_warm_phase=false;
	
        //=========================================================================
	// 4. Initialise flows

	outflow = outflow_types[params.parameters["flows"]["outflow"]["Form"]](params);

	// For inflow/radial flow we require rSFR
	double Rs = params.parameters["fundamentals"]["SolarRadius"];
	double ts = params.parameters["fundamentals"]["GalaxyAge"];
	double outflowrate_present = OutflowRate(Rs,ts,SFR(Rs,ts),GasReturnRate(Rs,ts));
	double prSFR=SFR(Rs,ts)-GasReturnRate(Rs,ts)+outflowrate_present;
	if(prSFR<0.)
		throw std::runtime_error("Reduced SFR<0 at R0 at final time\n");

	double K= params.parameters["fundamentals"]["Kennicutt-Schmidt_Coeff"];

	if (!check_param_given(params.parameters["fundamentals"],"Kennicutt-Schmidt_A",false)){
		double A = params.parameters["fundamentals"]["Kennicutt-Schmidt_A"];
		params.parameters["fundamentals"]["PresentGasDensitySun"]=pow(prSFR/A,1./K);
	}
	else{
		double PSG = params.parameters["fundamentals"]["PresentGasDensitySun"];
		params.parameters["fundamentals"]["Kennicutt-Schmidt_A"]=prSFR/pow(PSG,K);
	}


	inflow = inflow_types[params.parameters["flows"]["inflow"]["Form"]](
	            params,solar,prSFR);
	radialflow = radialflow_types[params.parameters["flows"]["radialflow"]["Form"]](params,solar,prSFR);
	if(gasdump)
		gasdumpflow = gasdump_types[params.parameters["flows"]["gasdump"]["Form"]](params,solar);
	else
		gasdumpflow = gasdump_types["None"](params,solar);
}
void Model::fill_initial_grids(void){

	double Zinit = params.parameters["fundamentals"]["InitialMetallicity"];
	double AlphaInit = extract_param(
	                    params.parameters["fundamentals"],"InitialAlpha",0.);

	Zinit *= solar->Z();
	metallicity->set_fixed_r_t_const(solar->Z());
	metallicity->set_fixed_t_const(Zinit,0);

	for(auto i=0u;i<mass_fraction.size();++i){
		auto initial_abundance = solar->scaled_solar_mass_frac(elements[i],
		                                                      Zinit,AlphaInit);
		mass_fraction[i].set_fixed_t_const(initial_abundance,0);
		if(use_warm_phase)
			mass_fraction_warm[i].set_fixed_t_const(initial_abundance,0);
	}

	// Fill gas and star initial grid
	VecDoub gas_radial_dist=gas_mass->grid_radial();
	VecDoub star_radial_dist=gas_mass->grid_radial();

	double K = params.parameters["fundamentals"]["Kennicutt-Schmidt_Coeff"];
	double A = params.parameters["fundamentals"]["Kennicutt-Schmidt_A"];
	double fracWarm=0.01;
	for(auto i=0u;i<gas_radial_dist.size();++i){
		star_radial_dist[i]=SFR(gas_mass->grid_radial()[i],0.);
		star_radial_dist[i]+=OutflowRate(gas_mass->grid_radial()[i],0.,
		                                 star_radial_dist[i],0.);
		gas_radial_dist[i]=pow(star_radial_dist[i]/A, 1./K);
	}
	if(warm_gas_mass){
		VecDoub warm_gas_radial_dist=gas_mass->grid_radial();
		for(auto i=0u;i<warm_gas_radial_dist.size();++i){
			warm_gas_radial_dist[i]=fracWarm*gas_radial_dist[i];
			if(warm_cooling_time>0.)
                            star_radial_dist[i]-=warm_gas_radial_dist[i]/warm_cooling_time;
			gas_radial_dist[i]=pow(star_radial_dist[i]/A, 1./K);
		}
		warm_gas_mass->set_fixed_t(warm_gas_radial_dist,0);
	}
	gas_mass->set_fixed_t(gas_radial_dist,0);
	reducedSFR->set_fixed_t(star_radial_dist,0);
	LOG(INFO)<<"Grids filled\n";
}

int Model::check_parameters(void){

	int err =0;
	// Check basics parameters
	auto F = params.parameters;
	std::vector<std::string> basic_blocks = {"fundamentals", "grids", "yields",
											 "typeIa","typeII", "flows",
											 "elements","data_folder"};
	for(auto s: basic_blocks) err+=check_param_given(F, s);
	if(err) return err;

	// Check fundamentals parameters
	std::vector<std::string> fund_blocks = {"MinimumMass",
											"MaximumMass",
											"GalaxyAge",
											"PresentSFR",
											"InitialMetallicity"};
	for(auto s: fund_blocks) err+=check_param_given(F["fundamentals"], s);
	if(err) return err;
	// Check grid parameters
	std::vector<std::string> grid_blocks = {"RadialGridPoints",
											"AgeGridPoints"};
	for(auto s: grid_blocks) err+=check_param_given(F["grids"], s);
	if(err) return err;

	if(params.parameters["grids"]["RadialGridPoints"]==1) single_zone=true;

	err+=check_param_given_matches_list(F["fundamentals"], imf_types, "IMF");
	err+=check_param_given_matches_list(F["fundamentals"], sfr_types, "SFR");
	err+=check_param_given_matches_list(F["fundamentals"], life_types,
	                                    "lifetimes");
	err+=check_param_given_matches_list_form(F, tia_types, "typeIa");
	err+=check_param_given_matches_list(F["yields"], typeII_types, "typeII");
	err+=check_param_given_matches_list(F["yields"], typeIa_types, "typeIa");
	err+=check_param_given_matches_list(F["yields"], agb_types, "AGB");
	// err+=check_param_given_matches_list(F["yields"], sagb_types, "SuperAGB");
	err+=check_param_given_matches_list_form(F["flows"], inflow_types,
	                                         "inflow");
	err+=check_param_given_matches_list_form(F["flows"], outflow_types,
	                                         "outflow");
	if(!single_zone){
		err+=check_param_given_matches_list_form(F, rm_types, "migration");
		err+=check_param_given_matches_list_form(F["flows"], radialflow_types,
		                                         "radialflow");
	}
	else{
		params.parameters["flows"]["radialflow"]["Form"]="None";
		params.parameters["migration"]["Form"]="None";
	}
	F = params.parameters["flows"];
	if (check_param_given(F, "gasdump", false)){
		gasdump=false;
		params.parameters["flows"]["gasdump"]["Form"]="None";
	}
	else{
		gasdump=true;
	        err+=check_param_given_matches_list_form(F, gasdump_types,
		                                         "gasdump");
	}
	F = params.parameters["elements"];
	for(auto f: F)
		if(element_index.find(f)==element_index.end()){
			LOG(INFO)<<"Element "<<f<<" not valid."<<std::endl;
		}

	if(params.parameters["grids"]["RadialGridPoints"]==1){
		single_zone=true;
		params.parameters["fundamentals"]["SolarRadius"]=1.;
		params.parameters["fundamentals"]["StarScaleLength"]=1.;
		params.parameters["grids"]["MinimumRadius"]=0.;
		params.parameters["grids"]["MaximumRadius"]=2.;
		params.parameters["migration"]["Form"]="None";
		params.parameters["flows"]["radialflow"]["Form"]="None";
		params.parameters["fundamentals"]["GasScaleLength"]=1;
	}
	else{
		err += check_param_given(params.parameters["grids"],"MinimumRadius");
		err += check_param_given(params.parameters["grids"],"MaximumRadius");
	}
	err+=check_param_given(params.parameters["fundamentals"],
	                       "Kennicutt-Schmidt_Coeff");
	if(check_param_given(params.parameters["fundamentals"],"PresentGasDensitySun",false)){
		err+=check_param_given(params.parameters["fundamentals"],
		                       "Kennicutt-Schmidt_A");
	}
	if(!check_param_given(params.parameters["fundamentals"],"Kennicutt-Schmidt_A",false) and
		!check_param_given(params.parameters["fundamentals"],"PresentGasDensitySun",false))
		LOG(INFO)<<"Both PresentGasDensitySun and Kennicutt-Schmidt_A given -- using Kennicutt-Schmidt_A";
	return err;
}

//=============================================================================
// Run models
int Model::step(unsigned nt, double dt){
	int err = 0;
	// Key Loop
	auto t = gas_mass->grid_time()[nt];
	auto NR = gas_mass->grid_radial().size();
	auto tp = t-.5*dt;

	auto gm_it=gas_mass->grid_fixed_t(nt);
	for(auto nR=0u;nR<NR;++nR) gm_it[nR]=1e90;
	std::vector<int> not_done(gm_it.size(),true);

	for(int N=0;N<iteratemax;++N){
        
        VecDoub migration_r(NR,0.);
	double migration_timestep=0.2;

	#pragma omp parallel for schedule(dynamic)
	for(auto nR=0u;nR<NR;++nR){
		// First we compute SFR and gas return rate -- this allows us
		// to find the gas needed in inflows in this time step
		auto R = gas_mass->grid_radial()[nR];
		// we set the metallicity to that of the previous time step so small
		// age interpolation works (i.e. doesn't use zero)
		metallicity->set(Z(R,t-dt),nR,nt);
		err+=check_metallicity(R, t, dt, nR, NR, nt);
                if(err) continue;
		auto starformrate = SFR(R,t);
		auto gas_return = GasReturnRate(R,t)*(1-warm_cold_ratio);
		auto outflowrate = OutflowRate(R,t,starformrate,gas_return);
		auto gasdumprate = GasDumpRate(R,t,1.);
	 	auto warmgasrate = 0.;
	        if(migration and nt>1)
                    migration_r[nR]=(rad_mig->convolve(gas_mass.get(),nR,nt,migration_timestep)-(*gas_mass)(nR,nt-1))/migration_timestep;
                if(use_warm_phase and warm_cooling_time>0.)
		    warmgasrate=(*warm_gas_mass)(nR,nt-1)/warm_cooling_time;
		auto rSFR = starformrate-gas_return
		            +outflowrate-warmgasrate-gasdumprate;
                rSFR+=-migration_r[nR];
                reducedSFR->set(rSFR,nR,nt);
	}
	if(err){return err;}
	#pragma omp parallel for schedule(dynamic)
	for(auto nR=0u;nR<NR;++nR){
		auto Xi=0.,R=0.;
		auto starformrate=0.,inflowrate=0.,outflowrate=0.,enrichrate=0.;
		double rad_flow_dm=0.,gas_return=0.,gas_dump_dm=0.;
		auto dmdt=0.,dmsdt=0.,sm_prev=0.,sm=0.,gm_prev=0.,gm=0.;
		auto rmlower=0.,rmupper=0.,rmlower_1=0.,rmupper_1=0.,gmprev_d=0.;
		auto area=0., area_up=0., area_down=0.,gmhere=0.,Xihere=0.;
		auto outer_gas_mass=0., e_outer_gas_mass=0.;
		auto warm_gm=0., warm_gm_prev =0.;

		R = gas_mass->grid_radial()[nR];
		
                // 1. Cold Gas
		gm_prev = (*gas_mass)(nR,nt-1);
		starformrate = SFR(R,tp);
		gas_return = GasReturnRate(R,tp);
		inflowrate = InflowRate(R,tp);
		outflowrate = OutflowRate(R,tp,starformrate,gas_return);
		rad_flow_dm = RadialFlowRateFromGrid(R, tp, dt, gm_prev,
		                                     nR, nt, NR, &err);
		// Gas dump evaluated at t (not tp)
		gas_dump_dm = GasDumpRate(R, t, dt);

		dmdt += -starformrate+gas_return*(1-warm_cold_ratio);
		dmdt += inflowrate;
		dmdt += -outflowrate+rad_flow_dm+gas_dump_dm;

		gm = gm_prev+dmdt*dt;

		if(gm<0.){
			starformrate=0.;rad_flow_dm=0.;outflowrate=0.;
			dmdt=gas_return+inflowrate-outflowrate;
                        gm=gm_prev+dmdt*dt;
		}

		// 2. Warm Gas
		if(use_warm_phase){
			warm_gm_prev = (*warm_gas_mass)(nR,nt-1);
			dmdt = gas_return*warm_cold_ratio + outflowrate*outflow_warm_fraction;
			warm_gm = warm_gm_prev + dmdt*dt;
                        if(warm_cooling_time>0.){
                            warm_gm += -warm_gm_prev*dt/warm_cooling_time;
			    gm += warm_gm_prev*dt/warm_cooling_time;
                        }
		}
                gm += migration_r[nR]*dt;
		gas_mass->set(gm,nR,nt);

		// 3. Stars
		sm_prev = (*stellar_mass)(nR,nt-1);
		dmsdt = starformrate;
		dmsdt -= gas_return;
		sm = sm_prev+dmsdt*dt;

		stellar_mass->set(sm,nR,nt);
		if(migration and nt>1){
			sm += rad_mig->convolve(stellar_mass.get(),nR,nt)-sm_prev;
			stellar_mass->set(sm,nR,nt);
		}

		if (use_warm_phase) warm_gas_mass->set(warm_gm, nR, nt);

		// 4. Enrichment
		for(auto e: elements){
			Xi = mass_fraction[e.first](nR,nt-1);
			dmdt = -starformrate*Xi;			  // star formation
			enrichrate=EnrichmentRate(e.second,R,tp);
			dmdt += enrichrate*(1-warm_cold_ratio); // re-enrich
			if(outflowrate>0.) dmdt -= OutflowRate(R,tp,starformrate*Xi,enrichrate);
			dmdt += inflowrate*inflow->Xi_inflow(e.second); // inflow
			if(rad_flow_dm!=0.) dmdt += EnrichRadialFlowRateFromGrid(e.second, R, tp, dt,
			                                                        gm_prev*Xi,
			                                                        nR, nt, NR, &err);
			dmdt += gas_dump_dm*gasdumpflow->Xi_inflow(e.second);
			double mass_now = Xi*gm_prev+dmdt*dt;
			double warm_mass_now = 0., Xiwarm=0.;

			if(use_warm_phase){
				Xiwarm = mass_fraction_warm[e.first](nR,nt-1);
				warm_mass_now = Xiwarm*warm_gm_prev;
				warm_mass_now += enrichrate*warm_cold_ratio*dt;
				// cooling
                                if(warm_cooling_time>0.){
                                    mass_now += Xiwarm*warm_gm_prev*dt/warm_cooling_time;
				    warm_mass_now -= Xiwarm*warm_gm_prev*dt/warm_cooling_time;
                                }
                                // inflow
				warm_mass_now += dt*inflowrate*(inflow->Xi_inflow(e.second)-Xiwarm);
                                mass_now += dt*inflowrate*(Xiwarm-inflow->Xi_inflow(e.second));
				// collecting outflow
				warm_mass_now += dt*outflow_warm_fraction*OutflowRate(R,tp,starformrate*Xi,enrichrate);

                  		mass_fraction_warm[e.first].set(warm_mass_now/warm_gm,nR,nt);
			}

			mass_fraction[e.first].set(mass_now/gm,nR,nt);

			if(migration and nt>1){
                               	mass_now+=(rad_mig->convolve_massfrac(gas_mass.get(),&mass_fraction[e.first],nR,nt,migration_timestep)-Xi*gm_prev)/migration_timestep*dt;
				mass_fraction[e.first].set(mass_now/gm,nR,nt);
			}
		}
		metallicity->set(1.-X(R,t)-Y(R,t),nR,nt);
		err = check_metallicity(R, t, dt, nR, NR, nt);

		if(fabs((gm-gm_it[nR])/gm)<tol) not_done[nR]=0;
		else not_done[nR]=1;
		}
	int done_sum=0;
	for_each(begin(not_done),end(not_done),
	         [&done_sum](int p){done_sum+=p;});
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
	auto NR = gas_mass->grid_radial().size();
	auto tp = t-.5*dt;

	for(auto nR=0u;nR<NR;++nR){
		double c = rad_mig->convolve(gas_mass.get(),nR,nt);
		gas_mass->set(c,nR,nt);
		std::cout<<gas_mass->grid_radial()[nR]<<" "<<t<<" "<<(*gas_mass)(nR,nt)<<std::endl;
	}
}

void Model::expand_grids(unsigned nt, double t){
	gas_mass->add_time(nt,t);
	stellar_mass->add_time(nt,t);
	reducedSFR->add_time(nt,t);
	metallicity->add_time(nt,t);
	for(unsigned mfn=0;mfn<mass_fraction.size();++mfn)
		mass_fraction[mfn].add_time(nt,t);
        if(use_warm_phase){
	    warm_gas_mass->add_time(nt,t);
	    for(unsigned mfn=0;mfn<mass_fraction_warm.size();++mfn)
		mass_fraction_warm[mfn].add_time(nt,t);
	}
}

void Model::run(void){
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
	metallicity->log10scale(solar->Z());
	for(auto e:elements)
		mass_fraction[e.first].log10scale(solar->mass_fraction(e.second));
}
int Model::check_metallicity(double R, double t, double dt, unsigned nR, unsigned NR, unsigned nt){
    if(Z(R,t)<0){
        metallicity->set(Z(R,t-dt),nR,nt);
        mass_fraction[element_index["H"]].set(mass_fraction[element_index["H"]](nR,nt-1),nR,nt);
        mass_fraction[element_index["He"]].set(mass_fraction[element_index["He"]](nR,nt-1),nR,nt);
        LOG(INFO)<<"Metallicity<0: trying to fix...: radial grid "<<nR<<" of "<<NR;
        LOG(INFO)<<"time "<<t<<std::endl;
        return 1;
    }
    else return 0;
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
	if(use_warm_phase)
		warm_gas_mass->write_hdf5(fout,"Mgas_warm");
	reducedSFR->write_hdf5(fout,"rSFR");
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
double Model::OutflowRate(double R, double t, double SFR, double GasReturn){return (*outflow)(R,t,SFR,GasReturn);}
double Model::GasDumpRate(double R, double t, double dt){
	return (*gasdumpflow)(R, t, dt);
}
// double Model::EnrichGasDumpRate(Element E, double R, double t, double dt){
// 	return gasdumpflow->elements(E, R, t, dt);
// }
double Model::RadialFlowRateFromGrid(double R, double t,
                                     double dt, double gm_prev,
                                     unsigned nR, unsigned nt,
                                     unsigned NR, int*err){
	if(single_zone)
		return 0.;
	auto outer_gas_mass=0.;
	auto Rdown=gas_mass->Rdown(nR);
	auto Rup=gas_mass->Rup(nR);

	if(nR<NR-1)
		outer_gas_mass=(*gas_mass)(nR+1,nt-1);
	else
		outer_gas_mass=gas_mass->log_extrapolate_high(Rup,nt-1);
	return RadialFlowRate(gm_prev,outer_gas_mass,R,Rdown,Rup,t,dt,err);
}
double Model::EnrichRadialFlowRateFromGrid(Element E, double R, double t,
                                           double dt, double e_gm_prev,
		                                   unsigned nR, unsigned nt,
		                                   unsigned NR, int*err){
	if(single_zone)
		return 0.;
	auto e_outer_gas_mass=0.;
	auto Rdown=gas_mass->Rdown(nR);
	auto Rup=gas_mass->Rup(nR);

	if(nR<NR-1)
		e_outer_gas_mass=(*gas_mass)(nR+1,nt-1)*mass_fraction[elements_r[E]](nR+1,nt-1);
	else
		e_outer_gas_mass=gas_mass->log_extrapolate_high(Rup,nt-1)*inflow->Xi_inflow(E);
	return RadialFlowRate(e_gm_prev,e_outer_gas_mass,R,Rdown,Rup,t,dt,err);
}

double Model::RadialFlowRate(double m, double mup, double R, double Rdown, double Rup, double t, double dt, int*err){
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
		double target_err=1e-3;
                return integrate(&_abundance_rate_2D,&P,target_err,0,"Divonne",&err,
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
	double integral = 0.;
	if(params.parameters["migration"]["Form"]=="None"){
		Rate_st P({this,R,t,which,mthresh});
		//GaussLegendreIntegrator GL(300);
		integrator GL(iterMAX_integral);
		integral = GL.integrate(&_gasreturn_rate,tmin,tmax,1e-4,&P);
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
		integral = integrate(&_gasreturn_rate_2D,&P,1e-3,0,"Divonne",&err,
						 "GasReturn "+which);
	}
	return integral;
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
