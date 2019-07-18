3d2
< unsigned iterMAX_integral=50000;
43a43,44
> 		if(e.empty())
> 			continue;
50c51,52
< 	// 4. Initialise flows
---
> 	// 4. Initialise flows -- we delay initialization of radial flows and
> 	//    inflow as they are sensitive to the initial metallicity choice.
52,56c54,56
< 	double Rs = params.parameters["fundamentals"]["SolarRadius"];
< 	double ts = params.parameters["fundamentals"]["GalaxyAge"];
< 	double prSFR=SFR(Rs,ts)-(1.-OutflowFraction(Rs,ts))*GasReturnRate(Rs,ts);
< 	inflow = inflow_types[params.parameters["flows"]["inflow"]["Form"]](params,prSFR);
< 	radialflow = radialflow_types[params.parameters["flows"]["radialflow"]["Form"]](params,prSFR);
---
> 	outer_metal =params.parameters["flows"]["radialflow"]["OuterEdgeMetallicity"];
> 	outer_edge =params.parameters["flows"]["radialflow"]["OuterEdge"];
> 	thresh = params.parameters["fundamentals"]["thresh"];
58a59
> 
62c63,65
<         metallicity->set_fixed_t_const(Zinit,0);
---
>     metallicity->set_const(Zinit);
> 
>     // Fill each element accounting for initial non-zero alpha enrichment
64c67
< 		auto initial_abundance = solar.scaled_solar_mass_frac(elements[i],Zinit);
---
> 		auto initial_abundance=solar.scaled_solar_mass_frac(elements[i],Zinit);
68a72
> 
76,84c80,85
< 	auto F = params.parameters["fundamentals"];
<         double A;
< 	if (F.find("Kennicutt-Schmidt_A") != F.end())
< 		A = F["Kennicutt-Schmidt_A"];
< 	else{
< 		double PresentGasDensity = params.parameters["fundamentals"]["PresentGasDensitySun"];
< 		double rSFR=SFR(Rs,ts)-(1.-OutflowFraction(Rs,ts))*GasReturnRate(Rs,ts);
< 		if(rSFR<0.)
< 			throw std::runtime_error("Reduced SFR<0 at R0 at final time\n");
---
> 	minAGBmass = params.parameters["yields"]["MinAGBMass"];
> 	double PresentGasDensity = params.parameters["fundamentals"]["PresentGasDensitySun"];
> 	double rSFR=SFR(Rs,ts)-(1.-OutflowFraction(Rs,ts))*GasReturnRate(Rs,ts);
> 	if(rSFR<0.)
> 		throw std::runtime_error("Reduced SFR<0 at R0 at final time\n");
> 	std::cout<<"Final reduced SFR = "<<rSFR<<std::endl;
86,87c87,93
< 		A=rSFR/pow(PresentGasDensity,KSN);
< 	}
---
> 	inflow = inflow_types[params.parameters["flows"]["inflow"]["Form"]](params,rSFR);
> 	radialflow = radialflow_types[params.parameters["flows"]["radialflow"]["Form"]](params,rSFR);
> 
> 	double A=rSFR/pow(PresentGasDensity,KSN);
> 	std::cout<<"Kennicutt-Schmidt coefficient A: "<<A<<std::endl;
> 
> 	double deltat = gas_mass->grid_time()[1]-gas_mass->grid_time()[0];
89c95,96
< 		star_radial_dist[i]=SFR(gas_mass->grid_radial()[i],0.);
---
> 		star_radial_dist[i]=SFR(gas_mass->grid_radial()[i],
> 		                        gas_mass->grid_time()[0]);
94,104d100
<         LOG(INFO)<<"Grids filled\n";
< 
<         if(gasdump){
< 		double mint=100000., tt;
<                 for(int i=0;i<gas_mass->grid_time().size();++i)
< 			if(fabs(gas_mass->grid_time()[i]-gasdumptime)<mint){
< 				tt = gas_mass->grid_time()[i];
< 				mint = fabs(gas_mass->grid_time()[i]-gasdumptime);
< 			}
< 		gasdumptime=tt;
< 	}
108,116d103
<         if(params.parameters["grids"]["RadialGridPoints"]==1){
<             single_zone=true;
<             params.parameters["fundamentals"]["SolarRadius"]=1.;
<             params.parameters["grids"]["MinimumRadius"]=0.;
<             params.parameters["grids"]["MaximumRadius"]=2.;
< 	    params.parameters["migration"]["Form"]="None";
< 	    params.parameters["flows"]["radialflow"]="None";
< 	    params.parameters["fundamentals"]["GasScaleLength"]=1;
<         }
158d144
<         if(!single_zone){
169d154
< 	}
182d166
<         if(!single_zone)
194,201d177
< 	F = params.parameters["flows"];
<         if (F.find("gasdump") != F.end()){
< 		gasdump=true;
< 		gasdumptime=F["gasdump"]["time"];
< 		gasdumpsurfacedensity=F["gasdump"]["surfacedensity"];
< 		gasdumpAlpha=F["gasdump"]["alpha"];
< 		gasdumpMetal=F["gasdump"]["metallicity"];
< 	}
225d200
< 		if(Z(R,t)<0.){ err=1; std::cerr<<"Returning here: "<<" "<<nR<<" "<<Z(gas_mass->grid_radial()[nR],t)<<" "<<Z(gas_mass->grid_radial()[nR],t-dt)<<std::endl;continue;}
229d203
< 		// if(nR==NR-1) reducedSFR->set(0.,nR,nt);
231,232c205
< 	if(err){return err;}
<         #pragma omp parallel for schedule(dynamic)
---
> 	#pragma omp parallel for schedule(dynamic)
239,240d211
< 		auto rmlower=0.,rmupper=0.,rmlower_1=0.,rmupper_1=0.,gmprev_d=0.;
< 		auto area=0., area_up=0., area_down=0.;
247,249d217
< 	        inflowrate = InflowRate(R,tp);
< 		//double A = params.parameters["fundamentals"]["Kennicutt-Schmidt_A"];
< 		//double K = params.parameters["fundamentals"]["Kennicutt-Schmidt_Coeff"];
251c219
< 		//inflowrate = (pow(starformrate/A,K)-gm_prev-(-starformrate+gas_return)*dt)/dt;
---
> 		inflowrate = InflowRate(R,tp);
253,254d220
< 		auto outer_gas_mass=0.;
< 		if(!single_zone){
257a224
> 		auto outer_gas_mass=0.;
260,264c227,231
< 		else
< 			outer_gas_mass=gas_mass->log_extrapolate_high(Rup,nt-1);
< 			// outer_gas_mass=exp(log((*gas_mass)(NR-1,nt-1))+(log((*gas_mass)(NR-1,nt-1))-log((*gas_mass)(NR-2,nt-1)))/(gas_mass->grid_radial()[NR-1]-gas_mass->grid_radial()[NR-2])*(Rup-gas_mass->grid_radial()[NR-1]));
< 		
< 		rad_flow_dm=RadialFlowRate(gm_prev,outer_gas_mass,R,Rdown,Rup,tp,dt,&err);
---
> 		else{
> 			if(outer_edge=="Extrap")
> 				outer_gas_mass=gas_mass->log_extrapolate_high(Rup,nt-1);
> 			else if (outer_edge=="Zero")
> 				outer_gas_mass=0.;
266c233,237
<                 dmdt += -starformrate;
---
> 		rad_flow_dm=RadialFlowRate(gm_prev,outer_gas_mass,R,Rdown,Rup,tp,dt,&err);
> 
> 		// rad_flow_dm=radialflow->radial_flow_rate(R,tp,reducedSFR.get());
> 
> 		dmdt += -starformrate;
270,276d240
< 		// if(nR==NR-1) dmdt=inflowrate+rad_flow_dm;
< 		if(((nR>NR-15) or (single_zone&&(nR==NR-1))) and fabs(gas_mass->grid_time()[nt]-gasdumptime)<1e-8 and gasdump){
< 			double mid,std=1.;
< 			if(single_zone) mid=gas_mass->grid_radial()[nR];
< 			else{mid=gas_mass->grid_radial()[NR-7]; std=gas_mass->grid_radial()[NR-5]-mid;}
< 			dmdt+=gasdumpsurfacedensity/dt*exp(-.5*pow((gas_mass->grid_radial()[nR]-mid)/std,2.));
< 		}
278,279c242,244
<                 if(gm<0.){starformrate=0.;rad_flow_dm=0.;dmdt=gas_return+inflowrate;gm=gm_prev+dmdt*dt;}
< 		
---
> 		if(fabs(dmdt*dt/gm_prev)>thresh)
> 			err=1;
> 
282,286d246
< 		/*if(gm<0.)
< 			throw std::runtime_error(
< 			           "More stars being formed than gas available in ring R="
< 			           +std::to_string(R));
< 		*/
288,313c248
< 		// if(migration){
< 		// 	area_up = gas_mass->annulus_area(nR+1);
< 		// 	area = gas_mass->annulus_area(nR);
< 		// 	area_down = gas_mass->annulus_area(nR-1);
< 		// 	rmlower=RadialMigrationKernel(Rdown,R,dt); // leaving gas moving in
< 		// 	rmupper=RadialMigrationKernel(Rup,R,dt); // leaving gas moving out
< 		// 	rmlower_1=RadialMigrationKernel(R,Rdown,dt)*area_down/area; // entering gas moving in
< 		// 	rmupper_1=RadialMigrationKernel(R,Rup,dt)*area_up/area; // entering gas moving out
< 		// 	if(nR==0)
< 		// 		gmprev_d=gas_mass->log_extrapolate_low(Rdown,nt-1);
< 		// 		// gmprev_d=exp(log((*gas_mass)(0,nt-1))+(log((*gas_mass)(1,nt-1))-log((*gas_mass)(0,nt-1)))/(gas_mass->grid_radial()[1]-gas_mass->grid_radial()[0])*(Rdown-gas_mass->grid_radial()[0]));
< 		// 	else
< 		// 		gmprev_d=(*gas_mass)(nR-1,nt-1);
< 		// 	// if(nR==0)
< 		// 	// 	gm+=outer_gas_mass*rmupper_1-gm_prev*rmupper;
< 		// 	// else
< 		// 		if(nR>0)
< 		// 		gm+=gmprev_d*rmlower_1+outer_gas_mass*rmupper_1-gm_prev*(rmupper+rmlower);
< 		// 	// std::cout<<R<<" "<<t<<" "<<rmlower<<" "<<rmupper<<" "<<rmlower_1<<" "<<rmupper_1<<" "<<" "<<Rdown<<" "<<Rup<<" "<<gm_prev<<" "<<gmprev_d<<" "<<gm_prev+dmdt*dt<<" "<<gm<<std::endl;
< 		// 		std::cout<<rmlower<<" "<<rmupper<<" "<<rmlower_1<<" "<<rmupper_1<<" "<<gmprev_d<<" "<<outer_gas_mass<<" "<<gm_prev<<" "<<Rdown<<" "<<Rup<<" "<<gas_mass->annulus_area(nR)<<" ";
< 		// }
< 
< 		std::cerr<<R<<" "<<t<<" "<<gm<<" "<<starformrate<<" "<<gas_return<<" "<<inflowrate<<" ";
< 		std::cerr<<rad_flow_dm<<" "<<reducedSFR->t_gradient(R,tp)<<" "<<outer_gas_mass<<" ";
< 		//std::cerr<<radialflow->gamma_g(R,Rdown,Rup,tp,dt,reducedSFR.get())<<" "<<radialflow->beta_g(R,Rdown,Rup,tp,dt,reducedSFR.get());
< 		//std::cerr<<" "<<radialflow->flow_rate((R+Rdown)*.5,tp,reducedSFR.get())<<" "<<radialflow->flow_rate((R+Rup)*.5,tp,reducedSFR.get());
---
> 		// std::cout<<R<<" "<<t<<" "<<gm<<" "<<starformrate<<" "<<gas_return<<" "<<inflowrate<<" "<<rad_flow_dm<<" "<<radialflow->radial_flow_rate(R,tp,reducedSFR.get())<<" "<<outer_gas_mass<<" "<<reducedSFR->t_gradient(R,tp)<<" "<<(*reducedSFR)(nR,nt)<<" "<<(*reducedSFR)(nR,nt-1)<<" "<<radialflow->flow_rate(R,tp,reducedSFR.get())<<" ";
323c258
< 		/*if(gm<0.)
---
> 		if(gm<0.)
327,328c262
< 		*/
< 		
---
> 
332a267,271
> 			if(outer_metal=="Fixed" and nR==NR-1){
> 				mass_fraction[e.first].set(mass_fraction[e.first](nR0,nt0)
> 				                           ,nR,nt);
> 				continue;
> 			}
339d277
< 			if(!single_zone){
343,344c281,286
< 			else
< 				e_outer_gas_mass *= mass_fraction[e.first](nR0,nt0);
---
> 			else{
> 				if(outer_metal=="Init")
> 					e_outer_gas_mass *= mass_fraction[e.first](nR0,nt0);
> 				else if(outer_metal=="Extrap")
> 					e_outer_gas_mass *= mass_fraction[e.first].log_extrapolate_high(Rup,nt-1);
> 			}
347,374d288
< 			}
< 
< 			if(((nR>NR-15) or (single_zone&&(nR==NR-1))) and fabs(gas_mass->grid_time()[nt]-gasdumptime)<1e-8 and gasdump){
<         			double metal = solar.Z()*pow(10.,gasdumpMetal);
<         			// Fill each element accounting for initial non-zero alpha enrichment
<                  		auto initial_abundance=solar.scaled_solar_mass_frac(e.second,metal);
<                                 if(is_alpha_element[e.second])
<                                     initial_abundance*=pow(10.,gasdumpAlpha);
<                         	double mid,std=1.;
< 				if(single_zone) mid=gas_mass->grid_radial()[nR];
< 				else{mid=gas_mass->grid_radial()[NR-7]; std=gas_mass->grid_radial()[NR-5]-mid;}
< 				dmdt+=gasdumpsurfacedensity*initial_abundance/dt*exp(-.5*pow((gas_mass->grid_radial()[nR]-mid)/std,2.));
<         		}
< 			// else dmdt=inflowrate*mass_fraction[e.first](nR0,nt0);
< 			// updated mass of element e
< 			// if(nR==NR-1) mass_fraction[e.first].set(Xi,nR,nt);
< 
< 			// double mig=0., lowerXi=0.;
< 			// if(migration){
< 			// 	if(nR==0)
< 			// 		lowerXi=mass_fraction[e.first](0,nt-1);
< 			// 	else
< 			// 		lowerXi=mass_fraction[e.first](nR-1,nt-1);
< 			// 	if(nR>0)
< 			// 	mig=gmprev_d*lowerXi*rmlower_1+e_outer_gas_mass*rmupper_1-gm_prev*(rmupper+rmlower)*Xi;
< 			// 	// if(nR==0)
< 			// 	// 	mig=e_outer_gas_mass*rmupper_1-gm_prev*rmupper*Xi;
< 			// }
376a291,292
> 			if(fabs(dmdt*dt/(Xi*gm))>thresh)
> 				err=1;
383,390c299
< 		if(Z(R,t)<0.) {
< 			metallicity->set(Z(R,t-dt),nR,nt); 
< 			mass_fraction[element_index["H"]].set(mass_fraction[element_index["H"]](nR,nt-1),nR,nt);
< 			mass_fraction[element_index["He"]].set(mass_fraction[element_index["He"]](nR,nt-1),nR,nt);
< 			err=1;std::cout<<"Problem "<<Z(R,t-dt)<<" "<<Z(R,t)<<std::endl;
< 		}
< 		// std::cout<<R<<" "<<t<<" "<<gm<<" "<<Z(R,t)<<std::endl;
< 		std::cerr<<" "<<X(R,t)<<" "<<Y(R,t)<<" "<<Z(R,t)<<std::endl;
---
> 		// std::cout<<TypeIIEnrichmentRate("He",R,tp)<<" "<<X(R,t)<<" "<<Y(R,t)<<" "<<Z(R,t)<<std::endl;
392,393c301
< 	if(err) return err;
< 	int iteratemax=1; auto gm_it=gas_mass->grid_fixed_t(nt);
---
> 	int iteratemax=3; auto gm_it=gas_mass->grid_fixed_t(nt);
414,415d321
< 			auto rmlower=0.,rmupper=0.,rmlower_1=0.,rmupper_1=0.,gmprev_d=0.;
< 			auto area=0., area_up=0., area_down=0.;
421a328
> 
423,424d329
< 		//double A = params.parameters["fundamentals"]["Kennicutt-Schmidt_A"];
< 		//double K = params.parameters["fundamentals"]["Kennicutt-Schmidt_Coeff"];
426,428d330
< 		//inflowrate = (pow(starformrate/A,K)-gm_prev-(-starformrate+gas_return)*dt)/dt;
< 			auto outer_gas_mass=0.;
< 			if(!single_zone){
431a334
> 			auto outer_gas_mass=0.;
433,435c336,342
< 				outer_gas_mass=(*gas_mass)(Rup,tp);
< 			else
< 				outer_gas_mass=gas_mass->log_extrapolate_high(Rup,tp);
---
> 				outer_gas_mass=(*gas_mass)(nR+1,nt-1);
> 			else{
> 				if(outer_edge=="Extrap")
> 					outer_gas_mass=gas_mass->log_extrapolate_high(Rup,nt-1);
> 				else if (outer_edge=="Zero")
> 					outer_gas_mass=0.;
> 			}
436a344
> 			gmhere = gm_prev;
439,440c347,349
< 			}
<                         dmdt += -starformrate;
---
> 			// rad_flow_dm=radialflow->radial_flow_rate(R,tp,reducedSFR.get());
> 
> 			dmdt += -starformrate;
444,450d352
< 
< 			if(((nR>NR-15) or (single_zone&&(nR==NR-1))) and fabs(gas_mass->grid_time()[nt]-gasdumptime)<1e-8 and gasdump){
<                         	double mid,std=1.;
< 				if(single_zone) mid=gas_mass->grid_radial()[nR];
< 				else{mid=gas_mass->grid_radial()[NR-7]; std=gas_mass->grid_radial()[NR-5]-mid;}
< 				dmdt+=gasdumpsurfacedensity/dt*exp(-.5*pow((gas_mass->grid_radial()[nR]-mid)/std,2.));
< 			}
452,478d353
<                         
<                         if(gm<0.){starformrate=0.;rad_flow_dm=0.;dmdt=gas_return+inflowrate;gm=gm_prev+dmdt*dt;}
< 
< 			// if(migration){
< 			// 	area_up = gas_mass->annulus_area(nR+1);
< 			// 	area = gas_mass->annulus_area(nR);
< 			// 	area_down = gas_mass->annulus_area(nR-1);
< 			// 	rmlower=RadialMigrationKernel(Rdown,R,dt); // leaving gas moving in
< 			// 	rmupper=RadialMigrationKernel(Rup,R,dt); // leaving gas moving out
< 			// 	rmlower_1=RadialMigrationKernel(R,Rdown,dt)*area_down/area; // entering gas moving in
< 			// 	rmupper_1=RadialMigrationKernel(R,Rup,dt)*area_up/area; // entering gas moving out
< 			// 	if(nR==0)
< 			// 		gmprev_d=gas_mass->log_extrapolate_low(Rdown,tp);
< 			// 	else
< 			// 		gmprev_d=(*gas_mass)(Rdown,tp);
< 			// 	// if(nR==0)
< 			// 	// 	gm+=outer_gas_mass*rmupper_1-gm_prev*rmupper;
< 			// 	// else
< 			// 	if(nR>0)
< 			// 		gm+=gmprev_d*rmlower_1+outer_gas_mass*rmupper_1-gmhere*(rmupper+rmlower);
< 			// 	std::cout<<rmlower<<" "<<rmupper<<" "<<rmlower_1<<" "<<rmupper_1<<" "<<gmprev_d<<" "<<outer_gas_mass<<" "<<gm_prev<<" "<<Rdown<<" "<<Rup<<" "<<gas_mass->annulus_area(nR)<<" ";
< 			// }
< 
< 			std::cerr<<R<<" "<<t<<" "<<gm<<" "<<starformrate<<" "<<gas_return<<" "<<inflowrate<<" ";
< 			std::cerr<<rad_flow_dm<<" "<<reducedSFR->t_gradient(R,tp)<<" "<<outer_gas_mass<<" ";
< 			//std::cerr<<radialflow->gamma_g(R,Rdown,Rup,tp,dt,reducedSFR.get())<<" "<<radialflow->beta_g(R,Rdown,Rup,tp,dt,reducedSFR.get())<<" ";
< 			//std::cerr<<radialflow->flow_rate((R+Rdown)*.5,tp,reducedSFR.get())<<" "<<radialflow->flow_rate((R+Rup)*.5,tp,reducedSFR.get());
487d361
< 				// std::cout<<gm-dmdt*dt-gm_prev<<" ";
491c365,366
< 			/*if(gm<0.)
---
> 			// std::cout<<R<<" "<<t<<" "<<gm<<" "<<starformrate<<" "<<gas_return<<" "<<inflowrate<<" "<<rad_flow_dm<<" "<<radialflow->radial_flow_rate(R,tp,reducedSFR.get())<<" "<<outer_gas_mass<<" "<<reducedSFR->t_gradient(R,tp)<<" "<<(*reducedSFR)(nR,nt)<<" "<<(*reducedSFR)(nR,nt-1)<<" "<<radialflow->flow_rate(R,tp,reducedSFR.get())<<" ";
> 			if(gm<0.)
495,496d369
< 			*/
< 		        //if(gm<0.) gm=0.;
508,509c381
<                          	if(!single_zone){
< 			        auto e_outer_gas_mass=outer_gas_mass;
---
> 				auto e_outer_gas_mass=outer_gas_mass;
511,513c383,389
< 					e_outer_gas_mass *= mass_fraction[e.first](Rup,tp);
< 				else
< 					e_outer_gas_mass *= mass_fraction[e.first](nR0,nt0);
---
> 					e_outer_gas_mass *= mass_fraction[e.first](nR+1,nt-1);
> 				else{
> 					if(outer_metal=="Init")
> 						e_outer_gas_mass *= mass_fraction[e.first](nR0,nt0);
> 					else if(outer_metal=="Extrap")
> 						e_outer_gas_mass *= mass_fraction[e.first].log_extrapolate_high(Rup,nt-1);
> 				}
516d391
< 				}
518,528d392
< 				if(((nR>NR-15) or (single_zone&&(nR==NR-1))) and fabs(gas_mass->grid_time()[nt]-gasdumptime)<1e-8 and gasdump){
<                                 	double metal = solar.Z()*pow(10.,gasdumpMetal);
<                                 	// Fill each element accounting for initial non-zero alpha enrichment
<                                 	auto initial_abundance=solar.scaled_solar_mass_frac(e.second,metal);
<                                 	if(is_alpha_element[e.second])
<                                 		initial_abundance*=pow(10.,gasdumpAlpha);
< 					double mid,std=1.;
< 					if(single_zone) mid=gas_mass->grid_radial()[nR];
< 					else{mid=gas_mass->grid_radial()[NR-7]; std=gas_mass->grid_radial()[NR-5]-mid;}
< 					dmdt+=gasdumpsurfacedensity*initial_abundance/dt*exp(-.5*pow((gas_mass->grid_radial()[nR]-mid)/std,2.));
< 				}
532d395
< 					// std::cout<<mass_now-dmdt*dt-Xi*gm_prev<<" ";
538,542c401
< 			if(Z(R,t)<0.) {
< 				metallicity->set(Z(R,tp),nR,nt); err=1;
< 			        mass_fraction[element_index["H"]].set(mass_fraction[element_index["H"]](nR,nt-1),nR,nt);
< 			        mass_fraction[element_index["He"]].set(mass_fraction[element_index["He"]](nR,nt-1),nR,nt);
< 			}
---
> 			// std::cout<<TypeIIEnrichmentRate("He",R,tp)<<" "<<X(R,t)<<" "<<Y(R,t)<<" "<<Z(R,t)<<std::endl;
544,545c403
< 			std::cerr<<" "<<X(R,t)<<" "<<Y(R,t)<<" "<<Z(R,t)<<std::endl;
< 			if(fabs((gm-gm_it[nR])/gm)<0.005) not_done[nR]=0;
---
> 			if(fabs((gm-gm_it[nR])/gm)<thresh/5.) not_done[nR]=0;
643,647c501
< 		int err;
< 		if(times[t]-times[t-1]>0.2)
< 			err=1;
< 		else
< 			err=step(t,times[t]-times[t-1]);
---
> 		int err=step(t,times[t]-times[t-1]);
673c527
< 	VecDoub sfr_grids,inflowss,snII_r,snIa_r;
---
> 	VecDoub sfr_grids,inflowss,snII_r,snIa_r,rsfr_grids;
675c529
<         for(auto t:gas_mass->grid_time()){
---
> 	for(auto t:gas_mass->grid_time()){
679a534
> 		rsfr_grids.push_back(rSFR(Rs,t));
684a540
> 	hdf5_write_1D_vector(fout,rsfr_grids,"rSFR");
702c558,561
< double Model::Z(double R, double t){ return (*metallicity)(R,t);}
---
> double Model::Z(double R, double t){
> 	// here we extrapolate for R<Rmin and R>Rmax
> 	return (*metallicity)(R,t,false,true);
> }
706a566
> double Model::rSFR(double R, double t){ return (*reducedSFR)(R,t);}
742a603
> // SNII rate using mass as variable -- below uses time
745c606,607
< 	return P->M->IMF(m)*P->M->SFR(P->R,P->t-P->M->Lifetime(m,0.02));
---
> 	return P->M->IMF(m)*P->M->SFR(P->R,
> 	                              P->t-P->M->Lifetime(m,P->M->Z(P->R,P->t)));
755c617
< 	Rate_st P(this,R,t,"",0.);
---
> 	Rate_st P(this,R,t,"",0.,0.);
766c628
<     if(P->which=="AGB" and M>P->mthresh)
---
>     if(P->which=="AGB" and (M>P->mthresh or M<P->minAGBmass))
779c641
<     if(tp>P->t) return 0.;
---
>     if(tp>P->t){ fval[0]=0.; return 0.;}
783,786c645,652
<     if(P->which=="AGB" and M>P->mthresh)
<     	return 0.;
<    	else if(P->which=="TypeII" and M<P->mthresh)
<    		return 0.;
---
>     if(P->which=="AGB" and (M>P->mthresh or M<P->minAGBmass)){
>     	fval[0]=0.;
>     	return 0;
>     }
>    	else if(P->which=="TypeII" and M<P->mthresh){
>     	fval[0]=0.;
>     	return 0;
>     }
795c661
< 		Rate_st P({this,R,t,which,mthresh});
---
> 		Rate_st P({this,R,t,which,mthresh,minAGBmass});
797c663
< 		integrator GL(iterMAX_integral);
---
> 		integrator GL(20000);
801,802c667,669
< 		double Rmin=0.,Rmax=30.;
< 		Rate_st_2D P(this,R,t,{tmin,Rmin},{tmax,Rmax},which,mthresh);
---
> 		double Rmin=0.,Rmax=params.parameters["grids"]["MaximumRadius"];
> 		Rmax+=10.;
> 		Rate_st_2D P(this,R,t,{tmin,Rmin},{tmax,Rmax},which,mthresh,minAGBmass);
818c685
< 	auto tmin = t-Lifetime(params.min_mass(),Z(R,0.));
---
> 	auto tmin = t-Lifetime(minAGBmass,Z(R,0.));
822c689
< 	if(tmax<0.) return 0.;
---
> 	if(tmax<0. or tmax<tmin) return 0.;
837c704
<     if(P->which=="AGB" and M>P->mthresh)
---
>     if(P->which=="AGB" and (M>P->mthresh or M<P->minAGBmass))
849c716
<     if(tp>P->t) return 0.;
---
>     if(tp>P->t){ fval[0]=0.; return 0.;}
853,856c720,727
<     if(P->which=="AGB" and M>P->mthresh)
<     	return 0.;
<    	else if(P->which=="TypeII" and M<P->mthresh)
<    		return 0.;
---
>     if(P->which=="AGB" and (M>P->mthresh or M<P->minAGBmass)){
>     	fval[0]=0.;
>     	return 0;
>     }
>    	else if(P->which=="TypeII" and M<P->mthresh){
>     	fval[0]=0.;
>     	return 0;
>     }
859,860d729
<     // if(P->which=="AGB")
<     // 	std::cout<<Rp<<" "<<tp<<" "<<P->E<<" "<<fval[0]<<std::endl;
864d732
< 
870c738
< 		EnrichRate_st P(this,R,t,E,which,mthresh);
---
> 		EnrichRate_st P(this,R,t,E,which,mthresh,minAGBmass);
873c741
< 		integrator GL(iterMAX_integral);
---
> 		integrator GL(20000);
877c745,755
< 		double Rmin=0.,Rmax=30.;
---
> 		double Rmin=0.,Rmax=params.parameters["grids"]["MaximumRadius"],Rmaxt,Rmint;
> 		Rmax+=10.;
> 		if(params.parameters["migration"]["Form"]=="Gaussian" or
> 		   params.parameters["migration"]["Form"]=="GaussianDrift"){
> 		   	double sR = (double)params.parameters["migration"]["sigmaR"];
> 		   	double tspread = 1.-tmin/(double)params.parameters["fundamentals"]["GalaxyAge"];
> 			Rmint=R-4.*sR*tspread;
> 			Rmin=(Rmint<Rmin?Rmin:Rmint);
> 			Rmaxt=R+4.*sR*tspread;
> 			Rmax=(Rmaxt>Rmax?Rmax:Rmaxt);
> 		}
880c758
< 		EnrichRate_st_2D P(this,R,t,E,xmin,xmax,which,mthresh);
---
> 		EnrichRate_st_2D P(this,R,t,E,xmin,xmax,which,mthresh,minAGBmass);
896c774
< 	auto tmin = t-Lifetime(params.min_mass(),Z(R,0.));
---
> 	auto tmin = t-Lifetime(minAGBmass,Z(R,0.));
900c778
< 	if(tmax<0.) return 0.;
---
> 	if(tmax<0. or tmax<tmin) return 0.;
920c798
<     if(P->which=="AGB" and M>P->mthresh)
---
>     if(P->which=="AGB" and (M>P->mthresh or M<P->minAGBmass))
932c810
<     if(tp>P->t) return 0.;
---
>     if(tp>P->t){ fval[0]=0.; return 0.;}
936,939c814,821
<     if(P->which=="AGB" and M>P->mthresh)
<     	return 0.;
<    	else if(P->which=="TypeII" and M<P->mthresh)
<    		return 0.;
---
>     if(P->which=="AGB" and (M>P->mthresh or M<P->minAGBmass)){
>     	fval[0]=0.;
>     	return 0;
>     }
>    	else if(P->which=="TypeII" and M<P->mthresh){
>     	fval[0]=0.;
>     	return 0;
>     }
941c823,824
<     fval[0]=P->M->TotalEjectedMass(M,Z)*P->M->IMF(M)*P->M->SFR(P->R,tp)*(-dmdt)*P->M->RadialMigrationKernel(P->R,Rp,age);
---
>     fval[0]=P->M->TotalEjectedMass(M,Z);
>     fval[0]*=P->M->IMF(M)*P->M->SFR(P->R,tp)*(-dmdt)*P->M->RadialMigrationKernel(P->R,Rp,age);
947,948c830,832
< 	if(params.parameters["migration"]["Form"]=="None"){
< 		Rate_st P({this,R,t,which,mthresh});
---
> 	double minAGBmassGR=minAGBmass;
> 	if(params.parameters["migration"]["Form"]=="None" or which=="TypeII"){
> 		Rate_st P({this,R,t,which,mthresh,minAGBmassGR});
950c834
< 		integrator GL(iterMAX_integral);
---
> 		integrator GL(20000);
954,960c838,857
< 		double Rmin=0.,Rmax=30.;
< 		// if(params.parameters["migration"]["Form"]=="Gaussian"){
< 		// 	Rmin=R-4.*(double)params.parameters["migration"]["sigmaR"];
< 		// 	Rmin=(Rmin<0.?0.:Rmin);
< 		// 	Rmax=R+4.*(double)params.parameters["migration"]["sigmaR"];
< 		// }
< 		Rate_st_2D P(this,R,t,{tmin,Rmin},{tmax,Rmax},which,mthresh);
---
> 		double Rmin=0.,Rmax=params.parameters["grids"]["MaximumRadius"],Rmaxt,Rmint;
> 		Rmax+=10.;
> 		if(params.parameters["migration"]["Form"]=="Gaussian" or
> 		   params.parameters["migration"]["Form"]=="GaussianDrift"){
> 		   	double sR = (double)params.parameters["migration"]["sigmaR"];
> 		   	double tspread = 1.-tmin/(double)params.parameters["fundamentals"]["GalaxyAge"];
> 			Rmint=R-4.*sR*tspread;
> 			Rmin=(Rmint<Rmin?Rmin:Rmint);
> 			Rmaxt=R+4.*sR*tspread;
> 			Rmax=(Rmaxt>Rmax?Rmax:Rmaxt);
> 		}
> 		// Rate_st_2D P(this,R,t,{tmin,R},{tmax,Rmax},which,mthresh,minAGBmassGR);
> 		// double err;
> 		// double integrand=integrate(&_gasreturn_rate_2D,&P,1e-3,0,"Divonne",&err,"GasReturn "+which);
> 		// P=Rate_st_2D(this,R,t,{tmin,Rmin},{tmax,R},which,mthresh,minAGBmassGR);
> 		// integrand+=integrate(&_gasreturn_rate_2D,&P,1e-3,0,"Divonne",&err,
> 		//                  "GasReturn "+which);
> 		// return integrand;
> 
> 		Rate_st_2D P(this,R,t,{tmin,Rmin},{tmax,Rmax},which,mthresh,minAGBmassGR);
962,963c859
< 		return integrate(&_gasreturn_rate_2D,&P,1e-3,0,"Divonne",&err,
< 		                 "GasReturn "+which);
---
> 		return integrate(&_gasreturn_rate_2D,&P,1e-3,0,"Divonne",&err,"GasReturn "+which);
972c868
< 	return SwitchGasReturnRate(R,t,tmin,tmax,"Type II");
---
> 	return SwitchGasReturnRate(R,t,tmin,tmax,"TypeII");
976c872
< 	auto tmin = t-Lifetime(params.min_mass(),Z(R,0.));
---
> 	auto tmin = t-Lifetime(minAGBmass,Z(R,0.));
993d888
< //=============================================================================
