#include "yields.h"
//=============================================================================
double yield_from_mass(unsigned e, double mej,double mrem,
                       double m, double Z, SolarAbundances &s){
	return (mej-(m-mrem)*s.scaled_solar_mass_frac(e,Z))/m;
}
double yield_from_mass(Element E, double mej,double mrem,
                       double m, double Z, SolarAbundances &s){
	return (mej-(m-mrem)*s.scaled_solar_mass_frac(E,Z))/m;
}
//=============================================================================
// Type Ia
//=============================================================================
TypeIaYields::TypeIaYields(ModelParameters M){
	// total ejected mass and remnant mass
	m_eject=Mch;
	m_remnant=0.;
}
Maeda_TypeIaYields::Maeda_TypeIaYields(ModelParameters M):TypeIaYields(M){
	std::string yields_folder=M.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("TypeIaYields file not found\n");
	}
	std::string line;
	for(int i=0;i<3;++i) std::getline(inFile,line);
	std::istringstream iss(line.substr(3));
	mass_eject=VecDoub((std::istream_iterator<double>(iss)),
	                    std::istream_iterator<double>());
}

Iwamoto_TypeIaYields::Iwamoto_TypeIaYields(ModelParameters M):TypeIaYields(M){
	std::string yields_folder=M.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("TypeIaYields file not found\n");
	}
	mass_eject = VecDoub(30,0.);
	std::string line;
	std::getline(inFile,line);
	std::getline(inFile,line);
	std::istringstream iss(line); std::string model;
	iss>>model;
	unsigned ncol=0;
	while(use_model.compare(model)!=0){
		iss>>model;ncol++;
	}
	VecDoub entries(6,0); std::string element;
	while(std::getline(inFile,line)){
		std::istringstream iss2(line);
		iss2>>element>>entries[0]>>entries[1]>>entries[2]
	      >>entries[3]>>entries[4]>>entries[5];
		element = element.substr(2);
		mass_eject[element_index[element]]+=entries[ncol-1];
	}
}
Seitenzahl_TypeIaYields::Seitenzahl_TypeIaYields(ModelParameters M):TypeIaYields(M){
	std::string yields_folder=M.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("TypeIaYields file not found\n");
	}
	mass_eject = VecDoub(32,0.);
	std::string line;
	std::getline(inFile,line);
	unsigned i=0;int atomic_number;std::string element;
	while(std::getline(inFile,line)){
		std::istringstream iss(line);
		iss>>element>>atomic_number>>mass_eject[i];
		i++;
	}
}
Thielemann_TypeIaYields::Thielemann_TypeIaYields(ModelParameters M)
:TypeIaYields(M){
	std::string yields_folder=M.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("TypeIaYields file not found\n");
	}
	mass_eject = VecDoub(41,0.);
	std::string line;
	for(auto i=0;i<3;++i)std::getline(inFile,line);
	std::istringstream iss(line);
	double tmp; iss>>tmp>>tmp;
	for(auto i=0;i<mass_eject.size();++i)
		iss>>mass_eject[i];
}

//=============================================================================
// General interface for Type II and AGB yields (metallicity dependent)
//=============================================================================
Yields_Single::Yields_Single(VecDoub mej,VecDoub ylds,double z, double m,double total_mass_eject, double mass_remnant):M(m),Z(z){
	mass_eject=mej; yields=ylds;
	m_remnant=mass_remnant; m_eject=total_mass_eject;
}
double Yields_Grid::interp_switch(Element E,int index,std::string s){
	if(s=="mass") return yields[index].mass(E);
	else if(s=="yield") return yields[index].yield(E);
	else if(s=="mass_remnant") return yields[index].mass_remnant();
	else if(s=="mass_ejected") return yields[index].mass_ejected();
	else throw std::invalid_argument(s+" is not a valid switch in interp");
}

double Yields_Grid::interp(Element E,double Mass, double Z,std::string m_or_y){

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
void Yields_Grid::sort_arrays(){
    std::vector<std::pair<unsigned,std::pair<double, double>>> met_mass;
    for(unsigned i=0;i<metal.size();++i)
        met_mass.push_back(std::pair<unsigned,std::pair<double, double>>({i,std::pair<double, double>({metal[i], mass_star[i]})}));
    sort(met_mass.begin(), met_mass.end(), comp_function_pair<double>);
    std::vector<size_t> order;for(unsigned i=0;i<metal.size();++i) order.push_back(met_mass[i].first);
    reorder(yields, order);
    sort(metal.begin(), metal.end());
    metal.erase( unique( metal.begin(), metal.end() ), metal.end() );
    sort(mass_star.begin(), mass_star.end());
    metal.erase( unique( metal.begin(), metal.end() ), metal.end() );
    mass_star.erase( unique( mass_star.begin(), mass_star.end() ), mass_star.end() );
}
//=============================================================================
// Type II
//=============================================================================
//=============================================================================
// Type II yields from Kobayashi 2006
TypeIIKobayashi2006::TypeIIKobayashi2006(ModelParameters M)
:TypeIIYields(M), solar(M){
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
		yields.push_back(Yields_Single(mass_eject,ylds,Z,m,finalm-mrem,mrem));
	}
    sort_arrays();
    // Now rescale
    rescale(&solar, solar_types[M.parameters["fundamentals"]["solar"]](M));
}
// Type II yields from Chieffi & Limongi 2004
TypeIIChieffiLimongi2004::TypeIIChieffiLimongi2004(ModelParameters M):TypeIIYields(M), solar(M){
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
		yields.push_back(Yields_Single(mass_eject,ylds,Z,m,meject,mrem));
	}
    sort_arrays();
    // Now rescale
    rescale(&solar, solar_types[M.parameters["fundamentals"]["solar"]](M));
}
// Type II yields from Chieffi & Limongi 2018
TypeIIChieffiLimongi2018::TypeIIChieffiLimongi2018(ModelParameters M):TypeIIYields(M), solar(M){
	std::string yields_folder=M.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("TypeIIChieffiLimongi2018 file not found\n");
	}
	std::string line; for(int i=0;i<51;++i) std::getline(inFile,line);
	int velp=0, fehp=0, vel, feh;
	std::vector<VecDoub> mass_eject(mass_list.size(),VecDoub(80,0.));
	VecDoub tmp(9,0.);
	std::string element;
	while(std::getline(inFile,line)){
		std::istringstream iss(line);
		iss>>vel>>feh>>element;
		element = std::regex_replace(element, std::regex(R"([\d])"), "");
		if(velp!=vel or fehp!=feh){
			if(velp==vel_target){
				for(int i=0;i<mass_eject.size();++i){
					double Z = metal_list[fehp];
					double meject = 0.;
					for(int j=0;j<mass_eject[i].size();++j)
						meject+=mass_eject[i][j];
					double m = mass_list[i];
					double mrem = m-meject;
					metal.push_back(log10(Z));
					mass_star.push_back(m);
					VecDoub ylds;
					for(int j=0;j<mass_eject[i].size();++j)
						ylds.push_back(yield_from_mass(i,mass_eject[i][j],
						                               mrem,
						                               m,Z,solar));
					yields.push_back(Yields_Single(mass_eject[i],ylds,
					                               Z,
					                               mass_star.back(),
					                               meject,
					                               mrem));
				}
				for(int i=0;i<mass_eject.size();++i)
					for(int j=0;j<mass_eject[i].size();++j)
						mass_eject[i][j]=0.;
			}
			velp=vel;fehp=feh;
		}
		for(int i=0;i<tmp.size();++i){
			iss>>tmp[i];
			mass_eject[i][element_index[element]]+=tmp[i];
		}
	}
    sort_arrays();
    // Now rescale
    rescale(&solar, solar_types[M.parameters["fundamentals"]["solar"]](M));
}
// Nugrid yields from Pignatari 2016 and Ritter 2018
// Also provide AGB yields
NugridYields::NugridYields(ModelParameters M):Yields_Grid(M), solar(M){
	std::string yields_folder=M.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("NuGrid yields file not found\n");
	}
	VecDoub mass_eject(83,0);
	std::string line;
	while(std::getline(inFile,line)){
		double Z, M, mrem;
		VecDoub mass_eject, ylds;
		if(line.substr(2,5).compare("Table")==0){
			int mindex = line.find("M"),
				 zindex = line.find("Z"),
				 endindex = line.find(")");
			M = std::stod(line.substr(mindex+2,zindex-1-(mindex+2)));
			Z = std::stod(line.substr(zindex+2,endindex-zindex-2));
			metal.push_back(log10(Z));mass_star.push_back(M);
			std::getline(inFile,line);std::getline(inFile,line);
			std::string element;
			std::istringstream iss(line);
			iss>>element>>element>>mrem;
			std::getline(inFile,line);
			double mass, minit, atomic_number;
			int i=0;
			while(std::getline(inFile,line)){
				std::istringstream iss2(line);
				iss2>>element>>mass>>minit>>atomic_number;
				mass_eject.push_back(mass);
				ylds.push_back(yield_from_mass(i,mass,mrem,
				                               M,Z,solar));
				i++;
				if(element.compare("Bi")==0)
					break;
			}
			yields.push_back(Yields_Single(mass_eject,ylds,Z,
			                                     M,M-mrem,mrem));
		}
	}
	sort_arrays();
    // Now rescale
    rescale(&solar, solar_types[M.parameters["fundamentals"]["solar"]](M));
}
//=============================================================================
// General interface to yields (metallicity dependent) where there isn't a
// regular grid of masses at each metallicity
//=============================================================================
double Yields_Grid_Irregular::interp_switch(Element E, unsigned iz,unsigned im, std::string s){
	if(s=="mass") return yields_vec[iz][im].mass(E);
	else if(s=="yield") return yields_vec[iz][im].yield(E);
	else if(s=="mass_remnant") return yields_vec[iz][im].mass_remnant();
	else if(s=="mass_ejected") return yields_vec[iz][im].mass_ejected();
	else throw std::invalid_argument(s+" is not a valid switch in interp");
}
double Yields_Grid_Irregular::interp(Element E,double Mass, double Z, std::string m_or_y){

	Z = log10(Z);
	int topmU,botmU,topmD,botmD,topz,botz;
	if(Z<metal.front()) Z=metal.front();
	if(Z>metal.back())  Z=metal.back();
    topbottom(metal,Z,&botz,&topz,"Z agb");

    auto mmx = mass_star_vec[topz].size()-1;
    if(Mass>mass_star_vec[topz][mmx]){ topmU=mmx; botmU=topmU-1;}
    else if(Mass<mass_star_vec[topz][0]){ topmU=1; botmU=0; }
    else topbottom(mass_star_vec[topz],Mass,&botmU,&topmU,"mass agb high z");

	auto yU1 = interp_switch(E,topz,botmU,m_or_y);
	auto yU2 = interp_switch(E,topz,topmU,m_or_y);

	mmx = mass_star_vec[botz].size()-1;
    if(Mass>mass_star_vec[botz][mmx]){ topmD=mmx; botmD=topmD-1;}
    else if(Mass<mass_star_vec[botz][0]){ topmD=1; botmD=0; }
    else topbottom(mass_star_vec[botz],Mass,&botmD,&topmD,"mass agb low z");

	auto yD1 = interp_switch(E,botz,botmD,m_or_y);
	auto yD2 = interp_switch(E,botz,topmD,m_or_y);

	if(yU1>0. and yU2>0. and yD1>0. and yD2>0.){
		yU1 = log10(yU1); yU2 = log10(yU2);
		yD1 = log10(yD1); yD2 = log10(yD2);
		auto yU=yU1+(Mass-mass_star_vec[topz][botmU])/(mass_star_vec[topz][topmU]-mass_star_vec[topz][botmU])*(yU2-yU1);
		auto yD=yD1+(Mass-mass_star_vec[botz][botmD])/(mass_star_vec[botz][topmD]-mass_star_vec[botz][botmD])*(yD2-yD1);
		return pow(10.,yD+(Z-metal[botz])/(metal[topz]-metal[botz])*(yU-yD));
	}
	else{
		auto yU=yU1+(Mass-mass_star_vec[topz][botmU])/(mass_star_vec[topz][topmU]-mass_star_vec[topz][botmU])*(yU2-yU1);
		auto yD=yD1+(Mass-mass_star_vec[botz][botmD])/(mass_star_vec[botz][topmD]-mass_star_vec[botz][botmD])*(yD2-yD1);
		return yD+(Z-metal[botz])/(metal[topz]-metal[botz])*(yU-yD);
	}
}
void Yields_Grid_Irregular::sort_arrays(void){

    std::vector<std::pair<unsigned,double>> met_;
    for(unsigned i=0;i<metal.size();++i)
        met_.push_back(std::pair<unsigned,double>({i,metal[i]}));
    sort(met_.begin(), met_.end(), comp_function<double>);
    std::vector<size_t> orderZ;
    for(unsigned i=0;i<metal.size();++i)
    	orderZ.push_back(met_[i].first);
    for(unsigned i=0;i<metal.size();++i){
    	std::vector<size_t> orderM;
    	std::vector<std::pair<unsigned,double>> mass_tmp;
    	for(unsigned j=0;j<mass_star_vec[i].size();++j)
	        mass_tmp.push_back(
	            std::pair<unsigned,double>({j, mass_star_vec[i][j]}));
	    sort(mass_tmp.begin(), mass_tmp.end(), comp_function<double>);
	    for(unsigned j=0;j<mass_star_vec[i].size();++j)
	    	orderM.push_back(mass_tmp[j].first);
	    reorder(mass_star_vec[i], orderM);
	    reorder(yields_vec[i], orderM);
    }
    reorder(yields_vec, orderZ);
    reorder(mass_star_vec, orderZ);
    reorder(metal, orderZ);
}
//=============================================================================
// AGB
//=============================================================================
// AGB yields from Karakas et al. (2010)
AGBYieldsKarakas::AGBYieldsKarakas(ModelParameters Params)
:Yields_Grid_Irregular(Params), solar(Params){
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
					yields_s[k]=yield_from_mass(k,mass_ejects[k],mrem,
					                            m,z,solar);
				yields_vec[iZ].push_back(Yields_Single(mass_ejects,yields_s,
				                                   z,m,m-mrem,mrem));
			}
			mass_ejects = VecDoub(30,0.);
			if(Z!=z){
				metal.push_back(log10(Z));
				mass_star_vec.push_back(VecDoub());
			    yields_vec.push_back(std::vector<Yields_Single>());
			    iZ+=1;
			}
			mass_star_vec[iZ].push_back(M);
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
	yields_vec[iZ].push_back(Yields_Single(mass_ejects,yields_s,Z,M,
	                                       M-mrem,mrem));

	sort_arrays();
    // Now rescale
    rescale(&solar,
            solar_types[Params.parameters["fundamentals"]["solar"]](Params));
}
//=============================================================================
// AGB yields from Ventura et al. (2013)
AGBYieldsVentura::AGBYieldsVentura(ModelParameters Params)
:Yields_Grid_Irregular(Params), solar(Params){
	std::string yields_folder=Params.parameters["data_folder"];
	std::ifstream inFile(yields_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+yields_folder+data_file<<std::endl;
		throw std::invalid_argument("AGBYieldsVentura file not found\n");
	}
	std::string line; std::getline(inFile,line); std::getline(inFile,line);
	VecDoub tmp(19,0.), tmp2(19,0.);
	std::map<unsigned, unsigned> index_element = {
		{0,0},{1,1},{2,2},{3,5},
		{4,5},{5,6},{6,7},{7,7},
		{8,7},{9,8},{10,9},{11,9},
		{12,10},{13,11},{14,11},
		{15,11},{16,12},{17,12},{18,13}};
	double mass, mrem, Z;
	while(std::getline(inFile,line)){
		if(line.substr(1,1)=="Z"){
			Z = std::stod(line.substr(3));
			metal.push_back(log10(Z));
			mass_star_vec.push_back(VecDoub());
		    yields_vec.push_back(std::vector<Yields_Single>());
		    std::getline(inFile,line);
		    std::getline(inFile,line);
			std::istringstream iss(line);
			for(auto i=0;i<tmp.size();++i)
				iss>>tmp[i];
			std::getline(inFile,line);
		}
		std::istringstream iss(line);
		iss>>mass;
		VecDoub mass_tmp(tmp.size(),0.), yields_tmp(tmp.size(),0.);
		for(auto i=0;i<tmp2.size();++i){
			iss>>tmp2[i];
			mass_tmp[index_element[i]]+=tmp2[i]+tmp[i]*mass;
		}
		iss>>mrem;
 		for(auto i=0;i<tmp2.size();++i){
			yields_tmp[index_element[i]]=yield_from_mass(i,
			                                        mass_tmp[index_element[i]],
			                                        mrem,mass,Z,solar);
		}
		mass_star_vec.back().push_back(mass);
		yields_vec.back().push_back(Yields_Single(mass_tmp,yields_tmp,
		                                          Z,mass,mass-mrem,mrem));
	}
	sort_arrays();
    // Now rescale
    rescale(&solar,
            solar_types[Params.parameters["fundamentals"]["solar"]](Params));
}

//=============================================================================
// Maps for creating yields instances from strings of class names
std::map<std::string, Yields_Grid*(*)(ModelParameters)> agb_types={
	{"Nugrid",&createAGBInstance<NugridYields>},
	{"Ventura",&createAGBInstance<AGBYieldsVentura>},
	{"Karakas",&createAGBInstance<AGBYieldsKarakas>},
	{"None",&createAGBInstance<Yields_Grid_None>}};

std::map<std::string, Yields_Grid*(*)(ModelParameters)> typeII_types={
	{"Nugrid",&createTypeIIInstance<NugridYields>},
	{"Kobayashi",&createTypeIIInstance<TypeIIKobayashi2006>},
	{"ChieffiLimongi2018",&createTypeIIInstance<TypeIIChieffiLimongi2018>},
	{"ChieffiLimongi2004",&createTypeIIInstance<TypeIIChieffiLimongi2004>},
	{"None",&createTypeIIInstance<Yields_Grid_None>}};

std::map<std::string, TypeIaYields*(*)(ModelParameters)> typeIa_types={
	{"Maeda",&createTypeIaInstance<Maeda_TypeIaYields>},
	{"Iwamoto",&createTypeIaInstance<Iwamoto_TypeIaYields>},
	{"Seitenzahl",&createTypeIaInstance<Seitenzahl_TypeIaYields>},
	{"Thielemann",&createTypeIaInstance<Thielemann_TypeIaYields>}};

//=============================================================================
// Set of yields
//     Combined AGB, Type Ia and Type II
//=============================================================================
YieldsSet::YieldsSet(ModelParameters M)
	:M(M),
	 agb(agb_types[M.parameters["yields"]["AGB"]](M)),
	 typeIa(typeIa_types[M.parameters["yields"]["typeIa"]](M)),
	 typeII(typeII_types[M.parameters["yields"]["typeII"]](M)){}

double YieldsSet::mass(Element E, double Mass, double Z){
	if(Mass<M.parameters["typeII"]["Min_typeII_SN_mass"])
		return agb->mass(E,Mass,Z);
	else
		return typeII->mass(E,Mass,Z);
}
double YieldsSet::mass_remnant(double Mass, double Z){
	if(Mass<M.parameters["typeII"]["Min_typeII_SN_mass"])
		return agb->mass_remnant(Mass,Z);
	else
		return typeII->mass_remnant(Mass,Z);
}
double YieldsSet::mass_ejected(double Mass, double Z){
	if(Mass<M.parameters["typeII"]["Min_typeII_SN_mass"])
		return agb->mass_ejected(Mass,Z);
	else
		return typeII->mass_ejected(Mass,Z);
}

double YieldsSet::typeIa_ejectedmass(Element E){
	return typeIa->mass(E);
}
//=============================================================================
