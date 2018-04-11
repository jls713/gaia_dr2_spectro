//=============================================================================
#include "iso_basti.h"

isochrone_johnson::isochrone_johnson(void){
    const unsigned N = 2000;
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

double isochrone_johnson::get_metallicity(std::vector<std::string> input_iso_l,
                                          std::string dir){
    auto input_iso = input_iso_l[0];
    int offset = dir.length()+1;
    // Load in isochrone file
    std::ifstream inFile;inFile.open(input_iso);
    if(!inFile.is_open()){
        std::cerr<<"Input isochrone file "<<input_iso
                 <<" won't open."<<std::endl;}
    // Record metallicity and age
    FeH = 0.06;
    if(input_iso.substr(2+offset,3).compare("sun") != 0){
        double Z_1 = atof(input_iso.substr(2+offset,1).c_str());
        double Z_2 = atof(input_iso.substr(3+offset,2).c_str());
        double Zin = Z_1*pow(10.,-Z_2);
        for(unsigned int i=0;i<BaSTI::FeHList.size();i++){
            FeH = BaSTI::FeHList[i];
            if(fabs(BaSTI::ZList[i]-Zin)<0.00001) break;
        }
    }
    inFile.close();
    return FeH;
}

void isochrone_johnson::fill(
			   std::vector<std::string> input_iso_list,
                           std::string dir, double Age){
    auto input_iso = input_iso_list[0];
    int offset = dir.length()+1;

    // offset gives the length of the string giving the isochrone directory

    // Load in isochrone file
    std::ifstream inFile;inFile.open(input_iso);
    if(!inFile.is_open()){std::cerr<<"Input isochrone file "<<input_iso<<" won't open."<<std::endl;}

    // Record metallicity and age
    FeH = 0.06;
    if(input_iso.substr(2+offset,3).compare("sun") != 0){
        double Z_1 = atof(input_iso.substr(2+offset,1).c_str());
        double Z_2 = atof(input_iso.substr(3+offset,2).c_str());
        double Zin = Z_1*pow(10.,-Z_2);
        for(unsigned int i=0;i<BaSTI::FeHList.size();i++){
            FeH = BaSTI::FeHList[i];
            if(fabs(BaSTI::ZList[i]-Zin)<0.00001) break;
        }
    }
    int ageindex = 13;
    if(FeH<-3.) ageindex+=2; // extra two 'ss' in filename
    age = atof(input_iso.substr(ageindex+offset,5).c_str())/1000.;

    std::string line;
    // Swallow header
    for(int i=0;i<8;i++)std::getline(inFile, line);

    // Load isochrone
    //(U-B)   (B-V)   (V-I)   (V-R)   (V-J)   (V-K)   (V-L)   (H-K)
    double IM,MASS,Mv,TEFF,LL,UB,BV,VI,VR,VJ,VK,VL,HK;
    int N=0;
    while(inFile >>IM>>MASS>>LL>>TEFF>>Mv>>UB>>BV>>VI>>VR
                >>VJ>>VK>>VL>>HK){
        InitialMass[N]=IM;Mass[N]=MASS;
        double J = Mv-VJ, K = Mv-VK, H = K+HK;
        VecDoub twoMASSColours = JKH_2MASS(J,K,H);
        L[N]=LL;
        mags["B"][N]=Mv+BV;
        mags["V"][N]=Mv;
        mags["Vp"][N]=VpMag(Mv+BV,Mv);
        mags["I"][N]=I2MASS(twoMASSColours[0],twoMASSColours[1]);
        mags["J"][N]=twoMASSColours[0];
        mags["H"][N]=twoMASSColours[2];
        mags["K"][N]=twoMASSColours[1];
        VecDoub JKHv = JKH_vista_from_2MASS(mags["J"][N],mags["K"][N],mags["H"][N]);
        mags["Jv"][N]=JKHv[0];
        mags["Hv"][N]=JKHv[2];
        mags["Kv"][N]=JKHv[1];
        Teff[N]=TEFF;
	Logg[N]=logg_calc(MASS,LL,TEFF);
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
    N_length=InitialMass.size();
    inFile.close();
}
