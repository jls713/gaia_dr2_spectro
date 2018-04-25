#include "extinction.h"

//=============================================================================
extinction_law::extinction_law(void){
	coeff_list["B"]=ABconst;
	coeff_list["V"]=AVconst;
	coeff_list["I"]=AIconst;
	coeff_list["J"]=AJconst;
	coeff_list["H"]=AHconst;
	coeff_list["K"]=AKconst;
}
double extinction_law::extinct_const(std::string band, double lTeff, double AV){
    if(band=="G")
	    return G_extinction(lTeff,AV);
	return coeff_list[band]*(1.-avgradient[band]*AV);
}
double extinction_law::av_gradient(std::string band){
	return avgradient[band];
}
double extinction_law::G_extinction(double lTeff, double AV){
	return AGconst*(1.-avgradient["G"]*AV);
}
double extinction_law::extinct_const_colour(std::vector<std::string> colour,
                                            double lTeff, double AV){
    return extinct_const(colour[0],lTeff,AV)-extinct_const(colour[1],lTeff,AV);
}
//=============================================================================
schlafly2017_extinction_law::schlafly2017_extinction_law(double RV){
	// Load extinction coeffs
	coeff_list.clear();
    std::string dir = parameters()["dir"]["extinction_coeffs"];
	std::ifstream inFile(dir+extinction_coeff_file);
	if(!inFile.good()){
		throw std::invalid_argument("coeff file not found\n");
	}
	VecDoub rvrange; std::string line; std::getline(inFile,line);
	std::istringstream iss(line);
	double rv;
	while(iss>>rv) rvrange.push_back(rv);
	std::map<std::string, VecDoub> coeffs;

	while(std::getline(inFile,line)){
		std::istringstream iss(line);
		VecDoub tmp; std::string mag;
		iss>>mag;
		double tmpd;
		while(iss>>tmpd) tmp.push_back(tmpd);
		coeffs[mag]=tmp;
		coeff_list[mag]=linterp(rvrange,coeffs[mag],RV,"constant");
	}
	inFile.close();

	// Transformed bands
	coeff_list["Vp"]=VpMag(coeff_list["B"],coeff_list["V"]);
	coeff_list["I"]=I2MASS(coeff_list["J"],coeff_list["K"]);
	VecDoub JKHv = JKH_vista_from_2MASS(coeff_list["J"],
	                                    coeff_list["K"],
	                                    coeff_list["H"]);
	coeff_list["Jv"]=JKHv[0];
	coeff_list["Hv"]=JKHv[2];
	coeff_list["Kv"]=JKHv[1];

	// Load G extinction

	std::ifstream inFile2(dir+G_extinction_coeff_file);
	if(!inFile2.good()){
		throw std::invalid_argument("G coeff file not found\n");
	}
	std::string line2; std::getline(inFile2,line2);
	std::istringstream iss2(line2);
	double teff;
	while(iss2>>teff) G_lTeff.push_back(teff);

	std::vector<VecDoub> Ggrid;

	while(std::getline(inFile2,line2)){
		std::istringstream iss3(line2);
		VecDoub tmp;
		double tmpd;
		while(iss3>>tmpd) tmp.push_back(tmpd);
		Ggrid.push_back(tmp);
	}
	for(double t: G_lTeff)
		G_extinct.push_back(linterp_2d(t,RV,G_lTeff,rvrange,Ggrid,"Gext"));
	inFile2.close();
	
        // Remove the last few elements of G_extinct as cause bad extrapolation
        G_extinct.pop_back();//G_extinct.pop_back();
        G_lTeff.pop_back();//G_lTeff.pop_back();
        // Normalize relative to A_V
	double VV = coeff_list["V"];
	std::cout<<VV<<std::endl;
        for (auto const& x : coeff_list){
		coeff_list[x.first]/=VV;
	}
	for (unsigned i=0; i<G_extinct.size(); ++i)
		G_extinct[i]/=VV;

    for(auto const &x: coeff_list)
        avgradient[x.first]=0.;
	std::ifstream inFile3(dir+avgradient_extinction_coeff_file);
	if(!inFile3.good()){
		throw std::invalid_argument("Av gradient file not found\n");
	}
    std::string band; double val;
    while(inFile3>>band>>val)
        avgradient[band]=val;

    inFile3.close();
    std::cout<<"Schlafly 2017 extinction law loaded\n";
}

double schlafly2017_extinction_law::G_extinction(double lTeff, double AV){
	// std::cout<<avgradient["G"]<<std::endl;
	return linterp(G_lTeff,G_extinct,lTeff,"linear")*(1.-avgradient["G"]*AV);
}
//=============================================================================
struct edensity_st{
	double l,b;
	sfd_extinction_map *EM;
	edensity_st(double l, double b, sfd_extinction_map *EM):l(l),b(b),EM(EM){}
};

static double edensity(double s, void*params){
	edensity_st* P = (edensity_st*)params;
	double SS = exp(s);
	return P->EM->density(P->l,P->b,SS)*SS;
}
//=============================================================================
double sfd_extinction_map::factor(double EBV){
	if(use_factor)
		return 0.6+0.2*(1-tanh((EBV-0.15)/0.3));
	else
		return 1.;
}

sfd_extinction_map::sfd_extinction_map(extinction_law *EL, VecDoub StS,
                               double RV, bool use_factor)
	:extinction_map(EL), StandardSolar(StS), RV(RV), use_factor(use_factor){
	Rfl = 1.12*StandardSolar[0];
	//=========================================================================
	// 1. Load in Healpix map
	Healpix_Ordering_Scheme scheme = NEST;
    HP_Map = Healpix_Map<double>();
    HP_Map.SetNside(NSIDE,scheme);

    read_Healpix_map_from_fits<double>(filename+"lambda_sfd_ebv.fits",HP_Map,1,2);
	//=========================================================================
	// 2. Create grids
    double deltal = 2.*PI/(double)(Nl);
    LL = create_range(deltal,2.*PI-deltal,Nl);
    double deltab = PI/(double)(Nb);
    BB = create_range(-PI/2.+deltab,PI/2.-deltab,Nb);
    SS = create_range(-3.6,2.5,Ns); //log s
	MAP=std::vector<std::vector<VecDoub>>(Nl,std::vector<VecDoub>(Nb,VecDoub(Ns,0.)));
	//=========================================================================
	// 3. Integrate along lines of sight
	#pragma openmp for parallel schedule(dynamic)
    for(unsigned i=0;i<LL.size();++i)
	    for(unsigned j=0;j<BB.size();++j){
	    	edensity_st EDD(LL[i],BB[j],this);
			double bs=PI/2.-BB[j];
			double EBV = HP_Map.interpolated_value(pointing(bs,LL[i]));
			double AVinf = RV*EBV*factor(EBV);
		    integrator INT(100000);
	    	double max = INT.integrate(&edensity,SS.front(),SS.back(),1e-3,&EDD);
	    	for(unsigned k=0;k<SS.size();++k){
		    	integrator INT2(100000);
		    	double now = INT2.integrate(&edensity,SS.front(),SS[k],1e-3,&EDD);
		    	MAP[i][j][k] = AVinf*now/max;
		    	if(MAP[i][j][k]<0.)
		    		std::cerr<<"A_V<0 for l="<<LL[i]<<", b="<<bs
		    				 <<", s="<<SS[k]<<std::endl;
	    	}
	    }

    std::cerr<<"Extinction map loaded"<<std::endl;
}
//=============================================================================
double sfd_extinction_map::evaluate(double l, double b, double s){

	double ls = log(s);
	if(l<0.)l+=2.*PI;
	int top_s,bot_s,top_l,bot_l,top_b,bot_b;
	//=========================================================================
	// 1. Find location in l
	if(l<LL.front() or l>LL.back()){ top_l=0; bot_l=LL.size()-1;}
	else topbottom<double>(LL,l,&bot_l,&top_l,"extinct l");
	double dl = (l-LL[bot_l])/(LL[top_l]-LL[bot_l]),mdl=1-dl;
	// Handle 2\pi periodicity of Galactic l properly
	if(l<LL.front()){
		dl = (l-LL[bot_l]+2.*PI)/(LL[top_l]-LL[bot_l]+2.*PI);mdl=1-dl;}
	else if(l>LL.back()){
		dl = (l-LL[bot_l])/(LL[top_l]-LL[bot_l]+2.*PI);mdl=1-dl;}

	//=========================================================================
	// 2. Find location in b
	if(b<BB.front()){bot_b=0;top_b=1;b=BB.front();}
	else if(b>BB.back()){top_b=BB.size()-1;bot_b=BB.size()-2;b = BB.back();}
	else topbottom<double>(BB,b,&bot_b,&top_b,"extinct b");
	double db = (b-BB[bot_b])/(BB[top_b]-BB[bot_b]),mdb=1-db;

	//=========================================================================
	// 3. Find location in s
	if(ls>SS.back())
		return dl*(db*MAP[top_l][top_b][Ns-1]+mdb*MAP[top_l][bot_b][Ns-1])+mdl*(db*MAP[bot_l][top_b][Ns-1]+mdb*MAP[bot_l][bot_b][Ns-1]);
	if(ls<SS.front())
		return dl*(db*MAP[top_l][top_b][0]+mdb*MAP[top_l][bot_b][0])+mdl*(db*MAP[bot_l][top_b][0]+mdb*MAP[bot_l][bot_b][0]);
	topbottom<double>(SS,ls,&bot_s,&top_s,"extinct s");
	double ds = (ls-SS[bot_s])/(SS[top_s]-SS[bot_s]),mds=1-ds;

	//=========================================================================
	// 4. 3D linear interpolation
	double result =
		   dl*(db*(ds*MAP[top_l][top_b][top_s]
	             +mds*MAP[top_l][top_b][bot_s])
			 +mdb*(ds*MAP[top_l][bot_b][top_s]
			     +mds*MAP[top_l][bot_b][bot_s]))
		 +mdl*(db*(ds*MAP[bot_l][top_b][top_s]
		         +mds*MAP[bot_l][top_b][bot_s])
		     +mdb*(ds*MAP[bot_l][bot_b][top_s]
		         +mds*MAP[bot_l][bot_b][bot_s]));
	if(result<0.){
		std::cerr<<"A_V<0 for l="<<l<<", b="<<b<<", s="<<s<<std::endl;
		std::cerr<<dl<<" "<<mdl<<" "<<db<<" "<<mdb<<" "<<ds<<" "<<mds<<std::endl;
	}
	return result;
}

double sfd_extinction_map::density(double l, double b, double s){

	VecDoub Pol = conv::GalacticToPolar({l,b,s},StandardSolar);
	double kfl = 1.+gammafl*MIN(Rfl,Pol[0]-Rfl);
	// Warp sinusoid in increasing phi coordinate in direction of Solar
	// motion -- which is minus our direction.
	double zw  = gammaw*MIN(Rw,Pol[0]-Rw)*sin(-Pol[1]);
	return exp((StandardSolar[0]-Pol[0])/hR-fabs(Pol[2]-zw)/kfl/hz);

}
//=============================================================================
combo_extinction_map::combo_extinction_map(extinction_law *EL)
	:extinction_map(EL){
	const std::string file = filename + "/combined_nside1024.fits";
    Healpix_Ordering_Scheme scheme = RING;//NEST;
    for(unsigned i=0;i<Nd;i++){
        HP_Map_grid.push_back(Healpix_Map<double>());
        HP_Map_grid[i].SetNside(NSIDE,scheme);
        read_Healpix_map_from_fits<double>(file,HP_Map_grid[i],i+1);
    }
    log10distance = create_range<double>(-1.2,1.8,Nd);
    std::cerr<<"Combo extinction map loaded"<<std::endl;
}

double combo_extinction_map::evaluate(double l, double b, double s,
                                      bool interp){
  double log10s = log10(s);
  std::cout<<l<<" "<<b<<std::endl;
  if(b>=PI/2. or b<=-PI/2.) return 0.;
  b=PI/2.-b;
  int bot_s, top_s;

  if(log10s<log10distance.front()){
    if(interp) return exp(HP_Map_grid[0].interpolated_value(pointing(b,l)));
    else return exp(HP_Map_grid[0][HP_Map_grid[0].ang2pix(pointing(b,l))]);
  }
  if(log10s>log10distance.back()){
    return exp(HP_Map_grid[Nd-1][HP_Map_grid[Nd-1].ang2pix(pointing(b,l))]);
  }
  topbottom<double>(log10distance,log10s,&bot_s,&top_s,"combo ext s");
  double yu,yd;
  if(interp){
    yu = HP_Map_grid[top_s].interpolated_value(pointing(b,l));
    yd = HP_Map_grid[bot_s].interpolated_value(pointing(b,l));
  }
  else{
    yu = HP_Map_grid[top_s][HP_Map_grid[top_s].ang2pix(pointing(b,l))];
    yd = HP_Map_grid[bot_s][HP_Map_grid[bot_s].ang2pix(pointing(b,l))];
  }
  return exp(yd+(log10s-log10distance[bot_s])/(log10distance[top_s]-log10distance[bot_s])*(yu-yd));
}
