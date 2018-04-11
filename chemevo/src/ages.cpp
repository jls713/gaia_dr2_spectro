#include "ages.h"
//=============================================================================
double MaederMeynet1989::operator()(double m, double Z){
	if(m<1.3) 			return pow(10.,-0.6545*log10(m)+1.);
	else if (m<=3.) 	return pow(10.,-3.7*log10(m)+1.35);
	else if (m<=7.) 	return pow(10.,-2.51*log10(m)+0.77);
	else if (m<=15.) 	return pow(10.,-1.78*log10(m)+0.17);
	else if (m<=60.) 	return pow(10.,-0.86*log10(m)-0.94);
	else 				return 1.2*pow(m,-1.85)+0.003;
}
double MaederMeynet1989::mass_star_dying_now(double t, double Z){
	if(t>8.42217) 		return pow(10.,-1.52788*(log10(t)-1.));
	else if(t>0.384283) return pow(10.,-0.27027*(log10(t)-1.35));
	if(t>0.0445455) 	return pow(10.,-0.3984*(log10(t)-.77));
	if(t>0.0119277) 	return pow(10.,-0.5618*(log10(t)-0.17));
	if(t>0.00339461) 	return pow(10.,-1.16279*(log10(t)+0.94));
	else 				return pow(0.83333*(t-0.003),-0.54054);
}
double MaederMeynet1989::dlogMdlogt(double t, double Z){
	if(t>8.42217) 		return -1.5279;
	else if(t>0.384283) return -0.27027;
	if(t>0.0445455) 	return -0.3984;
	if(t>0.0119277) 	return -0.5618;
	if(t>0.00339461) 	return -1.1628;
	else 				return -0.54054*t/(t-0.003);
}
//=============================================================================
// double Tinsley1980::operator()(double m, double Z){
// 	if(m<=1.) return 8.6;
// 	else
// }
double Kodama1997::operator()(double m, double Z){
	if(m<0.56)
		 return 50.;
	else if (m<=6.6)
		 return pow(10.,(0.334-sqrt(1.79-0.2232*(7.764-log10(m))))/0.1116);
	else return 1.2*pow(m,-1.85)+0.003;
}
double Kodama1997::mass_star_dying_now(double t, double Z){
	if (t>=0.04)
		 return pow(10.,7.764-4.48*(1.79-pow(0.334-0.1116*log10(t),2.)));
	else return pow(0.83333*(t-0.003),-0.54054);
}
double Kodama1997::dlogMdlogt(double t, double Z){
	if (t>=0.04) return -(0.334-0.1116*log10(t));
	else return -0.54054*t/(t-0.003);
}
//=============================================================================
double PadovaniMatteucci1993::operator()(double m, double Z){
	// Seems to be identical to Kodama but with t(m<0.56)=160Gyr!
	if(m<0.56)
		 return 160.;
	else if (m<=6.6)
		 return pow(10.,(0.334-sqrt(1.79-0.2232*(7.764-log10(m))))/0.1116);
	else return 1.2*pow(m,-1.85)+0.003;
}
double PadovaniMatteucci1993::mass_star_dying_now(double t, double Z){
	// Seems to be identical to Kodama but with t(m<0.56)=160Gyr!
	if (t>=0.04)
		 return pow(10.,7.764-4.48*(1.79-pow(0.334-0.1116*log10(t),2.)));
	else return pow(0.83333*(t-0.003),-0.54054);
}
double PadovaniMatteucci1993::dlogMdlogt(double t, double Z){
	if (t>=0.04) return -(0.334-0.1116*log10(t));
	else return -0.54054*t/(t-0.003);
}
//=============================================================================
struct find_mass_st{
	InterpolateAgeGrid *port;
	double t;
	double Z;
};
double find_mass(double m, void*p){
	find_mass_st *P=(find_mass_st *) p;
	return P->t-(*P->port)(pow(10.,m),P->Z);
}
Portinari1998::Portinari1998(ModelParameters M):
	InterpolateAgeGrid({-3.398,-2.398,-2.0969,-1.69897,-1.3010}){
	std::string data_folder = M.parameters["data_folder"];
	std::ifstream inFile(data_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+data_folder+data_file<<std::endl;
		throw std::invalid_argument("Portinari1998 data file not valid\n");
	}
	std::string line;
	std::getline(inFile,line);
	double m,tmp1,tmp2,tmp3,tmp4,tmp5;

	while(std::getline(inFile,line)){
		std::istringstream iss(line);
		iss>>m;mass_grid.push_back(log10(m));
		iss>>tmp1>>tmp2>>tmp3>>tmp4>>tmp5;
		lifetimes.push_back({log10(tmp1/1.e9),log10(tmp2/1.e9),
			log10(tmp3/1.e9),log10(tmp4/1.e9),log10(tmp5/1.e9)});
	}
	inFile.close();
	int2D = make_unique<interpolator2D>(mass_grid,metal_grid,lifetimes);

	root_find RF(1e-5,100000);
	auto max_age = lifetimes.front().back();
	auto min_age = lifetimes.back().back();
	age_grid = create_range(min_age,max_age,mass_grid.size()*4);
	unsigned na=0,nz=0;
	for(auto a: age_grid){
		masses.push_back(VecDoub(metal_grid.size(),0.));
		for(auto z:metal_grid){
			find_mass_st mm({this,pow(10.,a),pow(10.,z)});
			masses[na][nz]=RF.findroot(&find_mass,log10(0.3),log10(350.),&mm);
			++nz;
		}
		++na;nz=0;
	}
	int2D_back = std::make_shared<interpolator2D>(age_grid,metal_grid,masses);

}
double InterpolateAgeGrid::operator()(double m, double Z){
	Z = log10(Z);
	m = log10(m);

	if(Z<metal_grid.front())Z=metal_grid.front();
	if(Z>metal_grid.back()) Z=metal_grid.back();

	int topm,botm,topz,botz;
	if(m<mass_grid.front()){topm=1;botm=0;}
	else if(m>mass_grid.back()){botm=mass_grid.size()-2;topm=botm+1;}
    else return pow(10.,int2D->interpolate(m,Z));

    // Linear Extrapolation outside grid
    topbottom(metal_grid,Z,&botz,&topz);

	double ytt = lifetimes[topm][topz];
	double ytb = lifetimes[botm][topz];
	double ybt = lifetimes[topm][botz];
	double ybb = lifetimes[botm][botz];
	ybb = ybb+(Z-metal_grid[botz])/(metal_grid[topz]-metal_grid[botz])*(ytb-ybb);
	ybt = ybt+(Z-metal_grid[botz])/(metal_grid[topz]-metal_grid[botz])*(ytt-ybt);
	if(m>mass_grid.back()){
		// Extrapolation
		double t = ybt+(m-mass_grid[topm])/(mass_grid[topm]-mass_grid[botm])*(ybt-ybb);
		return pow(10.,t);
	}
	else if (m<mass_grid.front()){
		// Extrapolation
		double t = ybb+(m-mass_grid[botm])/(mass_grid[topm]-mass_grid[botm])*(ybt-ybb);
		return pow(10.,t);
	}
	else return pow(10.,ybb+(m-mass_grid[botm])/(mass_grid[topm]-mass_grid[botm])*(ybt-ybb));
}
double InterpolateAgeGrid::mass_star_dying_now(double t, double Z){
	Z = log10(Z);
	t = log10(t);
	int topm,botm,topz,botz;

	if(Z<metal_grid.front())Z=metal_grid.front();
	if(Z>metal_grid.back())Z=metal_grid.back();

	if(t<age_grid.front()){topm=1;botm=0;}
	else if(t>age_grid.back()){botm=age_grid.size()-2;topm=botm+1;}
    else return pow(10.,int2D_back->interpolate(t,Z));

    // Linear Extrapolation outside grid

    topbottom(metal_grid,Z,&botz,&topz);

	double ytt = masses[topm][topz];
	double ytb = masses[botm][topz];
	double ybt = masses[topm][botz];
	double ybb = masses[botm][botz];
	ybb = ybb+(Z-metal_grid[botz])/(metal_grid[topz]-metal_grid[botz])*(ytb-ybb);
	ybt = ybt+(Z-metal_grid[botz])/(metal_grid[topz]-metal_grid[botz])*(ytt-ybt);
	if(t>age_grid.back()){
		// Extrapolation
		double m = ybt+(t-age_grid[topm])/(age_grid[topm]-age_grid[botm])*(ybt-ybb);
		return pow(10.,m);
	}
	else if (t<age_grid.front()){
		// Extrapolation
		double m = ybb+(t-age_grid[botm])/(age_grid[topm]-age_grid[botm])*(ybt-ybb);
		return pow(10.,m);
	}
	else return pow(10.,ybb+(t-age_grid[botm])/(age_grid[topm]-age_grid[botm])*(ybt-ybb));
}

double InterpolateAgeGrid::dlogMdlogt(double t, double Z){
	Z = log10(Z);
	t = log10(t);
	int topm,botm,topz,botz;
	if(Z<metal_grid.front())Z=metal_grid.front();
	if(Z>metal_grid.back())Z=metal_grid.back();

	if(t<age_grid.front()){topm=1;botm=0;}
	else if(t>age_grid.back()){botm=age_grid.size()-2;topm=botm+1;}
    else return int2D_back->derivative_x(t,Z);

    // Linear Extrapolation outside grid
    topbottom(metal_grid,Z,&botz,&topz);

	double ytt = masses[topm][topz];
	double ytb = masses[botm][topz];
	double ybt = masses[topm][botz];
	double ybb = masses[botm][botz];

	ybb = (ytb-ybb)/(age_grid[topm]-age_grid[botm]);
	ybt = (ytt-ybt)/(age_grid[topm]-age_grid[botm]);

	return ybb+(Z-metal_grid[botz])/(metal_grid[topz]-metal_grid[botz])*(ybt-ybb);
}
//=============================================================================
PadovaIsochrones::PadovaIsochrones(ModelParameters M):
	InterpolateAgeGrid(create_range(-3., 0.7, 20)){
	std::string data_folder = M.parameters["data_folder"];
	std::ifstream inFile(data_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+data_folder+data_file<<std::endl;
		throw std::invalid_argument("PadovaIsochrones data file not valid\n");
	}
	std::string line;
	double m,tmp;

	while(std::getline(inFile,line)){
		std::istringstream iss(line);
		iss>>m;mass_grid.push_back(log10(m));
		VecDoub lifetime;
		while(iss>>tmp)
			lifetime.push_back(log10(tmp));
		lifetimes.push_back(lifetime);
	}
	inFile.close();
	int2D = make_unique<interpolator2D>(mass_grid,metal_grid,lifetimes);

	root_find RF(1e-5,100000);
	auto max_age = lifetimes.front().back();
	auto min_age = lifetimes.back().back();
	age_grid = create_range(min_age,max_age,mass_grid.size()*4);
	unsigned na=0,nz=0;
	for(auto a: age_grid){
		masses.push_back(VecDoub(metal_grid.size(),0.));
		for(auto z:metal_grid){
			find_mass_st mm({this,pow(10.,a),pow(10.,z)});
			masses[na][nz]=RF.findroot(&find_mass,log10(0.07),log10(350.),&mm);
			++nz;
		}
		++na;nz=0;
	}
	int2D_back = std::make_shared<interpolator2D>(age_grid,metal_grid,masses);

}
//=============================================================================
// Map for creating shared pointer instances of stellar lifetimes from string of class name
shared_map<StellarLifetime,ModelParameters> life_types ={
    {"MaederMeynet1989",
    	&createSharedInstance<StellarLifetime,MaederMeynet1989>},
    {"Kodama1997",
    	&createSharedInstance<StellarLifetime,Kodama1997>},
    {"PadovaniMatteucci1993",
    	&createSharedInstance<StellarLifetime,PadovaniMatteucci1993>},
    {"Portinari1998",
    	&createSharedInstance<StellarLifetime,Portinari1998>},
    {"PadovaIsochrones",
    	&createSharedInstance<StellarLifetime,PadovaIsochrones>}
};
//=============================================================================
