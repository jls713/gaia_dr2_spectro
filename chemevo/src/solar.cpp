#include "solar.h"
//=============================================================================
// These are atomic number - 1 upto Mo, and then we skip Tc so
// atomic number -2 beyond that, and then skip Pm so atomic number -3 beyond.
std::map<Element,int> element_index = {
	{"H",0},{"He",1},{"Li",2},{"Be",3},
	{"B",4},{"C",5},{"N",6},{"O",7},
	{"F",8},{"Ne",9},{"Na",10},{"Mg",11},
	{"Al",12},{"Si",13},{"P",14},{"S",15},
	{"Cl",16},{"Ar",17},{"K",18},{"Ca",19},
	{"Sc",20},{"Ti",21},{"V",22},{"Cr",23},
	{"Mn",24},{"Fe",25},{"Co",26},{"Ni",27},
	{"Cu",28},{"Zn",29},{"Ga",30},{"Ge",31},
	{"As",32},{"Se",33},{"Br",34},{"Kr",35},
	{"Rb",36},{"Sr",37},{"Y",38},{"Zr",39},
	{"Nb",40},{"Mo",41},{"Ru",42},{"Rh",43},
	{"Pd",44},{"Ag",45},{"Cd",46},{"In",47},
	{"Sn",48},{"Sb",49},{"Te",50},{"I",51},
	{"Xe",52},{"Cs",53},{"Ba",54},{"La",55},
	{"Ce",56},{"Pr",57},{"Nd",58},{"Sm",59},
	{"Eu",60},{"Gd",61},{"Tb",62},{"Dy",63},
	{"Ho",64},{"Er",65},{"Tm",66},{"Yb",67},
	{"Lu",68},{"Hf",69},{"Ta",70},{"W ",71},
	{"Re",72},{"Os",73},{"Ir",74},{"Pt",75},
	{"Au",76},{"Hg",77},{"Tl",78},{"Pb",79},
	{"Bi",80}
};
//=============================================================================
std::map<int,Element> element_index_reverse = {
	{0,"H"},{1,"He"},{2,"Li"},{3,"Be"},
	{4,"B"},{5,"C"},{6,"N"},{7,"O"},
	{8,"F"},{9,"Ne"},{10,"Na"},{11,"Mg"},
	{12,"Al"},{13,"Si"},{14,"P"},{15,"S"},
	{16,"Cl"},{17,"Ar"},{18,"K"},{19,"Ca"},
	{20,"Sc"},{21,"Ti"},{22,"V"},{23,"Cr"},
	{24,"Mn"},{25,"Fe"},{26,"Co"},{27,"Ni"},
	{28,"Cu"},{29,"Zn"},{30,"Ga"},{31,"Ge"},
	{32,"As"},{33,"Se"},{34,"Br"},{35,"Kr"},
	{36,"Rb"},{37,"Sr"},{38,"Y"},{39,"Zr"},
	{40,"Nb"},{41,"Mo"},{42,"Ru"},{43,"Rh"},
	{44,"Pd"},{45,"Ag"},{46,"Cd"},{47,"In"},
	{48,"Sn"},{49,"Sb"},{50,"Te"},{51,"I"},
	{52,"Xe"},{53,"Cs"},{54,"Ba"},{55,"La"},
	{56,"Ce"},{57,"Pr"},{68,"Nd"},{59,"Sm"},
	{60,"Eu"},{61,"Gd"},{62,"Tb"},{63,"Dy"},
	{64,"Ho"},{65,"Er"},{66,"Tm"},{67,"Yb"},
	{68,"Lu"},{69,"Hf"},{70,"Ta"},{71,"W "},
	{72,"Re"},{73,"Os"},{74,"Ir"},{75,"Pt"},
	{76,"Au"},{77,"Hg"},{78,"Tl"},{79,"Pb"},
	{80,"Bi"}
};
//=============================================================================
std::map<Element,bool> is_alpha_element = {
	{"H",false},{"He",false},{"Li",false},{"Be",false},
	{"B",false},{"C",true},{"N",false},{"O",true},
	{"F",false},{"Ne",true},{"Na",false},{"Mg",true},
	{"Al",false},{"Si",true},{"P",false},{"S",true},
	{"Cl",false},{"Ar",true},{"K",false},{"Ca",true},
	{"Sc",false},{"Ti",true},{"V",false},{"Cr",false},
	{"Mn",false},{"Fe",false},{"Co",false},{"Ni",false},
	{"Cu",false},{"Zn",false},{"Ga",false},{"Ge",false},
	{"As",false},{"Se",false},{"Br",false},{"Kr",false},
	{"Rb",false},{"Sr",false},{"Y",false},{"Zr",false},
	{"Nb",false},{"Mo",false},{"Ru",false},{"Rh",false},
	{"Pd",false},{"Ag",false},{"Cd",false},{"In",false},
	{"Sn",false},{"Sb",false},{"Te",false},{"I",false},
	{"Xe",false},{"Cs",false},{"Ba",false},{"La",false},
	{"Ce",false},{"Pr",false},{"Nd",false},{"Sm",false},
	{"Eu",false},{"Gd",false},{"Tb",false},{"Dy",false},
	{"Ho",false},{"Er",false},{"Tm",false},{"Yb",false},
	{"Lu",false},{"Hf",false},{"Ta",false},{"W ",false},
	{"Re",false},{"Os",false},{"Ir",false},{"Pt",false},
	{"Au",false},{"Hg",false},{"Tl",false},{"Pb",false},
	{"Bi",false}
};
//=============================================================================
AsplundSolarAbundances::AsplundSolarAbundances(ModelParameters M){
	std::string data_folder = M.parameters["data_folder"];
	std::ifstream inFile(data_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+data_folder+data_file<<std::endl;
		throw std::invalid_argument("Asplund data file not valid\n");
	}
	std::string line;
	for(auto i=0;i<4;++i) std::getline(inFile,line);
	double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;

	auto sumt=0.;
	while(std::getline(inFile,line)){
		std::istringstream iss(line);
		iss>>tmp1>>tmp2>>tmp3>>tmp4>>tmp5>>tmp6;
		masses.push_back(tmp2);
		abundances.push_back(tmp3);
		sumt+=masses.back()*pow(10.,abundances.back()-12.);
	}
	auto sumh=masses[0];
	auto sumhe=masses[1]*pow(10.,abundances[1]-12.);
	X_ = sumh/sumt; mfracs.push_back(X_);
	Y_ = sumhe/sumt; mfracs.push_back(Y_);
	Z_ = 1.-X_-Y_;

	for(auto i=2u;i<masses.size();++i)
		mfracs.push_back(masses[i]*pow(10.,abundances[i]-12.));
	inFile.close();
}
//=============================================================================
AndersTypeSolarAbundances::AndersTypeSolarAbundances(ModelParameters M, std::string data_file){
	std::string data_folder = M.parameters["data_folder"];
	std::ifstream inFile(data_folder+data_file);
	if(!inFile.good()){
		LOG(ERROR)<<"Could not open file "+data_folder+data_file<<std::endl;
		throw std::invalid_argument("Anders & Grevesse (1989) data file not valid\n");
	}
	std::string line;
	for(auto i=0;i<6;++i) std::getline(inFile,line);
	double tmp1,tmp2,tmp3;
	//char tmps = 'Mn';
	std::string tmps;
	auto sumt=0.;
	while(std::getline(inFile,line)){
		std::istringstream iss(line);
		iss>>tmp1>>tmps>>tmp2>>tmp3;
		masses.push_back(tmp3);
		abundances.push_back(tmp2);
		sumt+=masses.back()*pow(10.,abundances.back()-12.);
	}
	auto sumh=masses[0];
	auto sumhe=masses[1]*pow(10.,abundances[1]-12.);
	X_ = sumh/sumt; mfracs.push_back(X_);
	Y_ = sumhe/sumt; mfracs.push_back(Y_);
	Z_ = 1.-X_-Y_;

	for(auto i=2u;i<masses.size();++i)
		mfracs.push_back(masses[i]*pow(10.,abundances[i]-12.));
	inFile.close();
}
//=============================================================================
// Map for creating shared pointer instances of Solar from string of class name
shared_map< SolarAbundances,
            ModelParameters> solar_types ={
    {"Anders",&createSharedInstance<SolarAbundances,AndersSolarAbundances>},
    {"GrevesseNoel",&createSharedInstance<SolarAbundances,GrevesseNoelSolarAbundances>},
    {"GrevesseSauval",&createSharedInstance<SolarAbundances,GrevesseSauvalSolarAbundances>},
    {"Asplund",&createSharedInstance<SolarAbundances,AsplundSolarAbundances>}
};
//=============================================================================
