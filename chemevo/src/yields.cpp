#include "yields.h"
//=============================================================================
double remnant_mass(double initial_m){
	// Iben & Tutokov 1984b Pagel
	if(initial_m>9.5) return 1.5;
	else return 0.11*initial_m+0.45;
}
//=============================================================================
double yield_from_mass(unsigned e, double mej,double mrem,double m, double Z, SolarAbundances &s){
	return (mej-(m-mrem)*s.scaled_solar_mass_frac(e,Z))/m;
}
double yield_from_mass(Element E, double mej,double mrem,double m, double Z, SolarAbundances &s){
	return (mej-(m-mrem)*s.scaled_solar_mass_frac(E,Z))/m;
}
//=============================================================================
// Type Ia
//=============================================================================
TypeIaYields::TypeIaYields(ModelParameters M){
	std::string yields_folder=M.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("TypeIaYields file not found\n");
	}
	std::string line;
	for(int i=0;i<3;++i) std::getline(inFile,line);
	std::istringstream iss(line.substr(3));
	mass_eject=VecDoub((std::istream_iterator<double>(iss)),std::istream_iterator<double>());
	// total ejected mass and remnant mass
	m_eject=Mch;
	m_remnant=0.;
}

TypeIIYields_Single::TypeIIYields_Single(VecDoub mej,VecDoub ylds,double z, double m,double total_mass_eject, double mass_remnant):M(m),Z(z){
	mass_eject=mej; yields=ylds;
	m_remnant=mass_remnant; m_eject=total_mass_eject;
}
//=============================================================================
// Type II
//=============================================================================
double TypeIIYields::interp_switch(Element E,int index,std::string s){
	if(s=="mass") return yields[index].mass(E);
	else if(s=="yield") return yields[index].yield(E);
	else if(s=="mass_remnant") return yields[index].mass_remnant();
	else if(s=="mass_ejected") return yields[index].mass_ejected();
	else throw std::invalid_argument(s+" is not a valid switch in interp");
}

double TypeIIYields::interp(Element E,double Mass, double Z,std::string m_or_y){

	Z = log10(Z);
	int topm,botm,topz,botz;
	if(Z<metal.front()) Z=metal.front();
	if(Z>metal.back())  Z=metal.back();
    topbottom(metal,Z,&botz,&topz,"type II Z");

	if(Mass<mass_star.front()){botm=0;topm=1;}
	else if(Mass>mass_star.back()){botm=mass_star.size()-2;topm=botm+1;}
    else topbottom(mass_star,Mass,&botm,&topm,"type II mass");

	auto ytt = interp_switch(E,topz*mass_star.size()+topm,m_or_y);
	auto ytb = interp_switch(E,botz*mass_star.size()+topm,m_or_y);
	auto ybt = interp_switch(E,topz*mass_star.size()+botm,m_or_y);
	auto ybb = interp_switch(E,botz*mass_star.size()+botm,m_or_y);

	ytt=log10(ytt);ytb=log10(ytb);ybt=log10(ybt);ybb=log10(ybb);
	ybb = ybb+(Z-metal[botz])/(metal[topz]-metal[botz])*(ybt-ybb);
	ytb = ytb+(Z-metal[botz])/(metal[topz]-metal[botz])*(ytt-ytb);

	auto M = ybb+(Mass-mass_star[botm])/(mass_star[topm]-mass_star[botm])*(ytb-ybb);
	return pow(10.,M);
}
//=============================================================================
// Type II yields from Kobayashi 2006
TypeIIKobayashi2006::TypeIIKobayashi2006(ModelParameters M):TypeIIYields(M){
	std::string yields_folder=M.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("TypeIIKobayashi2006 file not found\n");
	}
	std::string line; std::getline(inFile,line);
	while(std::getline(inFile,line)){
		std::istringstream iss(line);
		double Z, m, mrem, finalm; iss>>m>>Z>>finalm>>mrem;
		if(Z==0.)Z+=small_metal;
		metal.push_back(log10(Z));mass_star.push_back(m);
		double mej;VecDoub mass_eject,ylds;
		auto i=0;
		while(iss>>mej){
			mass_eject.push_back(mej);
			ylds.push_back(yield_from_mass(i,mej,mrem,m,Z,solar));
			++i;
		}
		yields.push_back(TypeIIYields_Single(mass_eject,ylds,Z,m,finalm-mrem,mrem));
	}
	sort( metal.begin(), metal.end() );
	metal.erase( unique( metal.begin(), metal.end() ), metal.end() );
	sort( mass_star.begin(), mass_star.end() );
	mass_star.erase( unique( mass_star.begin(), mass_star.end() ), mass_star.end() );
}
// Type II yields from Chieffi & Limongi 2004
TypeIIChieffiLimongi2004::TypeIIChieffiLimongi2004(ModelParameters M):TypeIIYields(M){
	std::string yields_folder=M.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("TypeIIChieffiLimongi2004 file not found\n");
	}
	std::string line; std::getline(inFile,line);
	while(std::getline(inFile,line)){
		std::istringstream iss(line);
		double Z, m, mrem, meject; iss>>m>>Z>>meject>>mrem;
		if(Z==0.)Z+=small_metal;
		metal.push_back(log10(Z));mass_star.push_back(m);
		double mej;VecDoub mass_eject,ylds;
		auto i=0;
		while(iss>>mej){
			mass_eject.push_back(mej);
			ylds.push_back(yield_from_mass(i,mej,mrem,m,Z,solar));
			++i;
		}
		yields.push_back(TypeIIYields_Single(mass_eject,ylds,Z,m,meject,mrem));
	}
	sort( metal.begin(), metal.end() );
	metal.erase( unique( metal.begin(), metal.end() ), metal.end() );
	sort( mass_star.begin(), mass_star.end() );
	mass_star.erase( unique( mass_star.begin(), mass_star.end() ), mass_star.end() );
}

//=============================================================================
// Super AGB -- SUPER AGB NOT CURRENTLY IMPLEMENTED
//=============================================================================
// SuperAGBYields_Siess::SuperAGBYields_Siess(void){
// 	std::ifstream inFile(yields_folder+data_file);
// 	if(!inFile.good())
// 		std::cerr<<"Could not open file "+yields_folder+data_file<<std::endl;
// 	std::string line; for(int i=0;i<81;++i) std::getline(inFile,line);

// 	mass_star = {7.5,8.,8.5,9.,9.5,10.,10.5};

// 	while(std::getline(inFile,line)){
// 		std::istringstream iss(line);
// 		double Z, m; iss>>Z>>m;
// 		if(Z==0.)Z+=small_metal;
// 		metal.push_back(log10(Z));mass_star.push_back(m);
// 		yields.push_back(SuperAGBYields_Single(line));
// 	}
// }
//=============================================================================
// AGB
//=============================================================================
AGBYields_Single::AGBYields_Single(VecDoub mej,VecDoub ylds,double z, double m,double mass_remnant)
	:M(m),Z(z){
	mass_eject=mej; yields=ylds;
	m_remnant=mass_remnant; m_eject=m-m_remnant;
}

AGBYieldsKarakas::AGBYieldsKarakas(ModelParameters Params):AGBYields(Params){
	std::string yields_folder=Params.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("AGBYieldsKarakas file not found\n");
	}
	std::string line; for(int i=0;i<56;++i) std::getline(inFile,line);

	double z=0.,m=0.,mrem=0., Z=0., M,MREM;
	VecDoub yields_s(30,0.), mass_ejects(30,0.);
	int iZ=-1;

	while(std::getline(inFile,line)){
		std::istringstream iss(line);
		if(line.substr(6,1)=="y") continue;
		M = stod(line.substr(0,5));
		Z = stod(line.substr(9,6));
		MREM = stod(line.substr(17,5));
		if(M==m and Z==z and MREM!=mrem) continue;
		if(Z!=z or M!=m){
			if(m!=0. and z!=0.){
				for(unsigned k=0;k<yields_s.size();++k)
					yields_s[k]=yield_from_mass(k,mass_ejects[k],mrem,m,z,solar);
				yields[iZ].push_back(AGBYields_Single(mass_ejects,yields_s,z,m,mrem));
			}
			mass_ejects = VecDoub(30,0.);
			if(Z!=z){
				if(iZ==-1){
					metal.push_back(log10(Z));
					mass_star.push_back(VecDoub());
				    yields.push_back(std::vector<AGBYields_Single>());
				    iZ=0;
				}
				else{
					int botz, topz;
					if(log10(Z)<metal.front()) iZ=0;
					else if(log10(Z)>metal.back()) iZ=metal.size();
					else topbottom(metal,log10(Z),&botz,&topz);

					metal.insert(metal.begin()+iZ,log10(Z));
					mass_star.insert(mass_star.begin()+iZ,VecDoub());
				    yields.insert(yields.begin()+iZ,
				                  std::vector<AGBYields_Single>());
				}
			}
			mass_star[iZ].push_back(M);
		}
		m=M;z=Z;
		mrem = stod(line.substr(17,5));
		std::string el = line.substr(23,2);
		if(el.compare(1,1," ")==0)
			el=el[0];
		if(el=="g" or el=="n") continue;
		if(el=="p" or el=="D") el="H";
		mass_ejects[element_index[el]]+=stod(line.substr(42,8));
	}
	for(unsigned k=0;k<yields_s.size();++k)
		yields_s[k]=yield_from_mass(k,mass_ejects[k],mrem,m,z,solar);
	yields[iZ].push_back(AGBYields_Single(mass_ejects,yields_s,Z,M,mrem));
}
double AGBYieldsKarakas::interp_switch(Element E, unsigned iz,unsigned im, std::string s){
	if(s=="mass") return yields[iz][im].mass(E);
	else if(s=="yield") return yields[iz][im].yield(E);
	else if(s=="mass_remnant") return yields[iz][im].mass_remnant();
	else if(s=="mass_ejected") return yields[iz][im].mass_ejected();
	else throw std::invalid_argument(s+" is not a valid switch in interp");
}
double AGBYieldsKarakas::interp(Element E,double Mass, double Z, std::string m_or_y){

	Z = log10(Z);
	int topmU,botmU,topmD,botmD,topz,botz;
	if(Z<metal.front()) Z=metal.front();
	if(Z>metal.back())  Z=metal.back();
    topbottom(metal,Z,&botz,&topz,"Z agb");

    auto mmx = mass_star[topz].size()-1;
    if(Mass>mass_star[topz][mmx]){ topmU=mmx; botmU=topmU-1;}
    else if(Mass<mass_star[topz][0]){ topmU=1; botmU=0; }
    else topbottom(mass_star[topz],Mass,&botmU,&topmU,"mass agb high z");

	auto yU1 = interp_switch(E,topz,botmU,m_or_y);
	auto yU2 = interp_switch(E,topz,topmU,m_or_y);

	mmx = mass_star[botz].size()-1;
    if(Mass>mass_star[botz][mmx]){ topmD=mmx; botmD=topmD-1;}
    else if(Mass<mass_star[botz][0]){ topmD=1; botmD=0; }
    else topbottom(mass_star[botz],Mass,&botmD,&topmD,"mass agb low z");

	auto yD1 = interp_switch(E,botz,botmD,m_or_y);
	auto yD2 = interp_switch(E,botz,topmD,m_or_y);

	if(yU1>0. and yU2>0. and yD1>0. and yD2>0.){
		yU1 = log10(yU1); yU2 = log10(yU2);
		yD1 = log10(yD1); yD2 = log10(yD2);
		auto yU=yU1+(Mass-mass_star[topz][botmU])/(mass_star[topz][topmU]-mass_star[topz][botmU])*(yU2-yU1);
		auto yD=yD1+(Mass-mass_star[botz][botmD])/(mass_star[botz][topmD]-mass_star[botz][botmD])*(yD2-yD1);
		return pow(10.,yD+(Z-metal[botz])/(metal[topz]-metal[botz])*(yU-yD));
	}
	else{
		auto yU=yU1+(Mass-mass_star[topz][botmU])/(mass_star[topz][topmU]-mass_star[topz][botmU])*(yU2-yU1);
		auto yD=yD1+(Mass-mass_star[botz][botmD])/(mass_star[botz][topmD]-mass_star[botz][botmD])*(yD2-yD1);
		return yD+(Z-metal[botz])/(metal[topz]-metal[botz])*(yU-yD);
	}
}
//=============================================================================
// Maps for creating yields instances from strings of class names
std::map<std::string, AGBYields*(*)(ModelParameters)> agb_types={
	{"Karakas",&createAGBInstance<AGBYieldsKarakas>},
	{"None",&createAGBInstance<AGBYields>}};
std::map<std::string, SuperAGBYields*(*)(ModelParameters)> sagb_types={
	{"None",&createSAGBInstance<SuperAGBYields_None>}
	// ,{"Siess",&createSAGBInstance<SuperAGBYields_Siess>}
	};
std::map<std::string, TypeIIYields*(*)(ModelParameters)> typeII_types={
	{"Kobayashi",&createTypeIIInstance<TypeIIKobayashi2006>},
	{"ChieffiLimongi",&createTypeIIInstance<TypeIIChieffiLimongi2004>}};
std::map<std::string, TypeIaYields*(*)(ModelParameters)> typeIa_types={
	{"Maeda",&createTypeIaInstance<TypeIaYields>}};
//=============================================================================
// Set of yields
//=============================================================================
YieldsSet::YieldsSet(ModelParameters M)
	:M(M),
	 agb(agb_types[M.parameters["yields"]["AGB"]](M)),
	 super_agb(sagb_types[M.parameters["yields"]["SuperAGB"]](M)),
	 typeIa(typeIa_types[M.parameters["yields"]["typeIa"]](M)),
	 typeII(typeII_types[M.parameters["yields"]["typeII"]](M)){}

double YieldsSet::mass(Element E, double Mass, double Z){
	if(typeid(*super_agb)==typeid(SuperAGBYields_None)){
		if(Mass<M.parameters["typeII"]["Min_typeII_SN_mass"])
			return agb->mass(E,Mass,Z);
		else
			return typeII->mass(E,Mass,Z);
	}
	throw std::invalid_argument("Mass not possible using Super AGB "+std::string(typeid(*super_agb).name()));
    return 0.;
}
double YieldsSet::mass_remnant(double Mass, double Z){
	if(typeid(*super_agb)==typeid(SuperAGBYields_None)){
		if(Mass<M.parameters["typeII"]["Min_typeII_SN_mass"])
			return agb->mass_remnant(Mass,Z);
		else
			return typeII->mass_remnant(Mass,Z);
	}
	throw std::invalid_argument("Mass remnant not possible using Super AGB "+std::string(typeid(*super_agb).name()));
    return 0.;
}
double YieldsSet::mass_ejected(double Mass, double Z){
	if(typeid(*super_agb)==typeid(SuperAGBYields_None)){
		if(Mass<M.parameters["typeII"]["Min_typeII_SN_mass"])
			return agb->mass_ejected(Mass,Z);
		else
			return typeII->mass_ejected(Mass,Z);
	}
	throw std::invalid_argument("Mass ejected not possible using Super AGB "+std::string(typeid(*super_agb).name()));
    return 0.;
}

double YieldsSet::typeIa_ejectedmass(Element E){
	return typeIa->mass(E);
}
//=============================================================================
