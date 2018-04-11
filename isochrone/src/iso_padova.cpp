//=============================================================================
#include "iso_padova.h"
//=============================================================================

void find_feh_age(std::ifstream &inFile,std::string input_iso,double Age,double *age, double *FeH){
    bool found=0;
    std::string line;
    const std::string age_st = "Age =";
    int count=0;
    while (std::getline(inFile, line)){
    count++;
    if (line.find(age_st) != std::string::npos){
            *age = pow(10.,Age);
        if(fabs(atof(line.substr(line.find(age_st)+7,10).c_str())-*age)/(*age)<0.01){
            found=1;
                *age /=1e9;
            *FeH = atof(line.substr(45,6).c_str());
                break;
            }
        }
    }
    if(found==0){std::cerr<<"Search failed -- Padova isochrone"<<" "<<Age<<" "<<input_iso<<std::endl; return;}
    std::getline(inFile, line);
}
// //=============================================================================

isochrone_padova::isochrone_padova(void){
    const unsigned N = 500;
    InitialMass = VecDoub(N,0.);
    Mass = VecDoub(N,0.);
    L = VecDoub(N,0.);
    Teff = VecDoub(N,0.);
    Logg = VecDoub(N,0.);
    es = VecDoub(N,0.);
    std::vector<std::string> maglist =
        {"B","V",
         "u","g","r","i","z",
         "gP","rP","iP","zP",
         "G","GBP","GRP",
         "J","H","K",
         "W1","W2","W3","W4",
         "Jv","Hv","Kv","Vp","I"};
    for(auto m: maglist)
    mags[m]=VecDoub(N,0.);
    mag_maps ={
    {"g",{std::make_pair("G",8),
          std::make_pair("GBP",9),std::make_pair("GRP",10)}},
    {"s",{std::make_pair("u",8),std::make_pair("g",9),std::make_pair("r",10),
          std::make_pair("i",11),std::make_pair("z",12)}},
    {"u",{std::make_pair("B",9),std::make_pair("V",10)}},
    {"2",{std::make_pair("J",8),std::make_pair("H",9),
          std::make_pair("K",10),
          std::make_pair("W1",18),std::make_pair("W2",19),
          std::make_pair("W3",20),std::make_pair("W4",21)}},
    {"p",{std::make_pair("gP",8),
          std::make_pair("rP",9),std::make_pair("iP",10),
          std::make_pair("zP",11)}}};
}

double isochrone_padova::get_metallicity(std::vector<std::string> input_iso,
                                       std::string dir){
    double f =0.;
    std::ifstream inFile;inFile.open(input_iso[0]);
    std::string line;
    const std::string age_st = "Age =";
    while (std::getline(inFile, line)){
        if (line.find(age_st) != std::string::npos){
            f = atof(line.substr(45,6).c_str());
            break;
        }
    }
    inFile.close();
    return f;
}

void isochrone_padova::fill(std::vector<std::string> input_iso,
                            std::string dir, double Age){
    bool first=true;
    for(auto s: input_iso){
    std::ifstream inFile;inFile.open(s);
    if(!inFile.is_open()){
        std::cerr<<"Input isochrone file "<<s<<" won't open."<<std::endl;
    }
    find_feh_age(inFile,s,log10(Age)+9.,&age,&FeH);
    VecDoub X;std::string line;
    unsigned N=0;
    while (std::getline(inFile, line)){
        if (line.find("#") != std::string::npos or line.empty()) break;
        X = string2vec<double>(line);
        if(first){
            InitialMass[N]=X[2];
            Mass[N]=X[3];
            L[N]=X[4];
            Teff[N]=X[5];
            Logg[N]=X[6];
        }
        for(auto n:mag_maps[s.substr(dir.length()+1,1)])
            mags[n.first][N]=X[n.second];
        N++;
    }
    N_length=N;
    inFile.close();
    if(first){
        InitialMass.resize(N);
        Mass.resize(N);
        L.resize(N);
        Teff.resize(N);
        Logg.resize(N);
        es.resize(N);
        first=false;
    }
    for(auto n:mag_maps[s.substr(dir.length()+1,1)]){
        if(mags[n.first].size()<N)
            std::cout<<"Mags longer:"<<n.first<<" "<<s<<std::endl;
        mags[n.first].resize(N);
    }
    double mM = maxmass();mf = InitialMass;
    for(unsigned i=0;i<N_length;++i)
        mf[i]=mf[i]/mM;
    }
    for(unsigned N=0;N<N_length;++N){
        mags["Vp"][N] = VpMag(mags["B"][N],mags["V"][N]);
        mags["I"][N] = I2MASS(mags["J"][N],mags["K"][N]);
        VecDoub JKHv = JKH_vista_from_2MASS(mags["J"][N],
                                            mags["K"][N],
                                            mags["H"][N]);
        mags["Jv"][N]=JKHv[0];
        mags["Hv"][N]=JKHv[2];
        mags["Kv"][N]=JKHv[1];
    }
}
