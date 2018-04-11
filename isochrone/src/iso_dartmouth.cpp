//=============================================================================
#include "iso_dartmouth.h"
//=============================================================================

isochrone_dartmouth::isochrone_dartmouth(void){
    const unsigned N = 1000;
    InitialMass = VecDoub(N,0.);
    Mass = VecDoub(N,0.);
    L = VecDoub(N,0.);
    Teff = VecDoub(N,0.);
    Logg = VecDoub(N,0.);
    es = VecDoub(N,0.);
    std::vector<std::string> maglist =
        {"B","V",
         "J","H","K",
         "Jv","Hv","Kv","Vp","I"};
    for(auto m: maglist)
    mags[m]=VecDoub(N,0.);
}
double isochrone_dartmouth::get_metallicity(std::vector<std::string> input_iso_list, std::string dir){
    double ff = 0.;
    auto input_iso = input_iso_list[0];
    // Load in isochrone file
    std::ifstream inFile;inFile.open(input_iso);
    if(!inFile.is_open()){
        std::cerr<<"Input isochrone file "<<input_iso<<" won't open."<<std::endl;
    }
    std::string line;
    while (std::getline(inFile, line))
    {
        if (line.find("Fe/H") != std::string::npos){
            std::getline(inFile, line);
            VecDoub X = string2vec<double>(line.erase(0,1));
            ff = X[4];
            break;
        }
    }
    return ff;
}

void isochrone_dartmouth::fill(std::vector<std::string> input_iso_list,
                               std::string dir, double Age){
    auto input_iso = input_iso_list[0];
    // Load in isochrone file
    std::ifstream inFile;inFile.open(input_iso);
    if(!inFile.is_open()){
        std::cerr<<"Input isochrone file "<<input_iso<<" won't open."<<std::endl;
    }
    std::string line; bool found=0;
    const std::string age_st = "#AGE=";
    while (std::getline(inFile, line) and found==0)
    {
        if (line.find("Fe/H") != std::string::npos){
            std::getline(inFile, line);
            VecDoub X = string2vec<double>(line.erase(0,1));
            FeH = X[4];
        }
	if (line.find(age_st) != std::string::npos){
            if(atof(line.substr(5,6).c_str())==Age){
                found=1;break;
            }
        }
    }
    if(found==0) std::cerr<<"Search failed -- Dartmouth isochrone"<<" "<<Age<<" "<<input_iso<<std::endl;
    std::getline(inFile, line);
    VecDoub X;int N=0;
    while (std::getline(inFile, line)){
        if (line.find("#") != std::string::npos or line.empty()) break;
        X = string2vec<double>(line);
        InitialMass[N]=X[1];
        Mass[N]=X[1];
        L[N]=X[4];
        Teff[N]=X[2];
        Logg[N]=X[3];

        mags["B"][N]=X[6];
        mags["V"][N]=X[7];
        mags["Vp"][N]=VpMag(X[6],X[7]);
        mags["J"][N]=X[10];
        mags["H"][N]=X[11];
        mags["K"][N]=X[12];
        mags["I"][N]=I2MASS(mags["J"][N],mags["K"][N]);
        VecDoub JKHv = JKH_vista_from_2MASS(mags["J"][N],mags["K"][N],mags["H"][N]);
        mags["Jv"][N]=JKHv[0];
        mags["Hv"][N]=JKHv[2];
        mags["Kv"][N]=JKHv[1];
        ++N;
    }
    InitialMass.resize(N);
    Mass.resize(N);
    L.resize(N);
    Teff.resize(N);
    Logg.resize(N);
    es.resize(N);
    for(auto n: mags)
        mags[n.first].resize(N);
    age = Age;
    N_length=InitialMass.size();
    inFile.close();
}
//=============================================================================
