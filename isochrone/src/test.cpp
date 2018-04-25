//=============================================================================
// Test programs
//=============================================================================
#include "distances.h"
#include <map>
//=============================================================================
#include "gtest/gtest.h"
const double test_err = 1e-7;

//=============================================================================

namespace {
//=============================================================================

TEST(prior,prior){
    double l = 20.*PI/180.;
    double b = 8.*PI/180.;
    new_prior_2018 GP(conv::StandardSolarPAUL);
    for(auto r=0.1;r<20.;r+=0.1)
    std::cout<<r<<" "
             <<GP.prior({8.2-r*cos(l)*cos(b),-r*sin(l)*cos(b),r*sin(b)},0.,10.)
             <<std::endl;
}

template<class isochrone_g>
void test_distance(DistanceCalculator<isochrone_g> *D,
                   isochrone_grid<isochrone_g> *iso,
                   galaxy_prior *GP,
                   extinction_law *EL,
                   VecInt LOC_AFM={45,10,751}){

    int Na = LOC_AFM[0], NF = LOC_AFM[1], Nm = LOC_AFM[2];
    VecDoub lb = {0.4,1.};
    double s = 10.;
    double Teff = iso->iso(NF,Na)->logTeff(Nm);
    double Z = iso->iso(NF,Na)->feh();
    double age = iso->iso(NF,Na)->tau();
    double mass = iso->iso(NF,Na)->initial_mass(Nm);
    double logg = iso->iso(NF,Na)->logg(Nm);
    double logg_err = 0.1, feh_err = 0.2, teff_err=0.05;
    VecDoub covar = {feh_err*feh_err,0.,0.,
                     teff_err*teff_err,0.,
                     logg_err*logg_err};

    // return 0;
    double DM = 5.*log10(100.*s), AV=0.;
    std::vector<std::string> maglist = {"J","H","K","G","GRP","GBP","Jv"};
    VecDoub mag, magerr;

    for(auto s: maglist){
        magerr.push_back(0.022);
        mag.push_back(DM+iso->iso(NF,Na)->mag(Nm,s)+EL->extinct_const(s)*AV);
    }

    std::cout<<age<<" "<<mass<<" ";
    printVector(maglist);
    printVector(mag);
    std::cout<<Z<<" "<<Teff<<" "<<logg<<std::endl;

    printVector(D->prob_distance(mag,
                                Z,Teff,logg,
                                lb[0],lb[1],
                                magerr,
                                covar,
                                GP,
                                false, // don't return full pdf
                                5.,
                                maglist));
    double Aprior=0., parallax=0., parallax_error=-1.;

    sfd_extinction_map EM(EL,conv::StandardSolarPAUL);
    printVector(D->prob_distance_extinct(mag,
                                Z,Teff,logg,
                                lb[0],lb[1],
                                magerr,
                                covar,
                                GP,
                                false, // don't return full pdf
                                5.,
                                maglist,
                                &EM,
                                Aprior,
                                parallax,
                                parallax_error));
}

TEST(BaSTI,johnson){
    isochrone_grid<isochrone_johnson> iso("BaSTI",1,0.5);
    DistanceCalculator<isochrone_johnson> D(&iso);
    new_prior_2018 GP(conv::StandardSolarPAUL);
    schlafly2017_extinction_law EL;
    test_distance(&D,&iso,&GP,&EL);
}
TEST(Padova,all){
    isochrone_grid<isochrone_padova> iso("Padova",1,0.5);
    DistanceCalculator<isochrone_padova> D(&iso);
    new_prior_2018 GP(conv::StandardSolarPAUL);
    schlafly2017_extinction_law EL;
    test_distance(&D,&iso,&GP,&EL,{30,6,145});
}
TEST(Dartmouth,Dartmouth){
    isochrone_grid<isochrone_dartmouth> iso("Dartmouth");
    DistanceCalculator<isochrone_dartmouth> D(&iso);
    new_prior_2018 GP(conv::StandardSolarPAUL);
    schlafly2017_extinction_law EL;
    test_distance(&D,&iso,&GP,&EL,{5,4,15});
}
TEST(ExtMap,ExtMap){
    schlafly2017_extinction_law EL;
    sfd_extinction_map EM(&EL,conv::StandardSolarPAUL);
}
TEST(MaxMass,BaSTI){
    isochrone_grid<isochrone_johnson> iso("BaSTI");
    for(double a=1.;a<9.;a+=0.02)
        std::cout<<a<<" "<<iso.max_mass(0.,a)<<" "<<iso.max_mass(-0.1,a)<<" "<<iso.max_mass(0.1,a)<<std::endl;
}
TEST(MaxAge,BaSTI){
    isochrone_grid<isochrone_johnson> iso("BaSTI");
    for(double M=2.;M>0.8;M-=0.01)
        std::cout<<M<<" "<<iso.max_age(0.,M)<<" "<<iso.max_age(-0.1,M)<<" "<<iso.max_age(0.1,M)<<std::endl;
}
TEST(ES,BaSTI){
    isochrone_grid<isochrone_johnson> iso("BaSTI");
    int NA=34;
    double age = iso.iso(9,NA)->tau();
    double Z = iso.iso(9,NA)->feh();
    double age1 = iso.iso(9,NA+1)->tau();
    double Z1 = iso.iso(9,NA+1)->feh();
    for(unsigned i=0;i<iso.iso(9,NA)->N();++i){
        double mass = iso.iso(9,NA)->initial_mass(i);
        double tmax = iso.max_age(Z,iso.iso(9,NA)->initial_mass(i));
        std::cout<<mass<<" "<<age/tmax<<" "<<iso.iso(9,NA)->logTeff(i)<<" "<<iso.iso(9,NA)->logg(i)<<" ";
        mass = iso.iso(9,NA+1)->initial_mass(i);
        tmax = iso.max_age(Z1,iso.iso(9,NA+1)->initial_mass(i));
        std::cout<<mass<<" "<<age1/tmax<<" "<<iso.iso(9,NA+1)->logTeff(i)<<" "<<iso.iso(9,NA+1)->logg(i)<<" ";
        mass = iso.iso(9,NA)->initial_mass(i);
        tmax = iso.max_age(Z1,iso.iso(9,NA+1)->initial_mass(i));
        std::cout<<mass<<" "<<(age+0.25)/tmax<<" ";
        VecDoub ff = iso.interp_es(Z1,age+0.25,mass,{"V"});
        std::cout<<ff[0]<<" "<<ff[1]<<std::endl;
    }
}
TEST(MaxMass,Padova){
    isochrone_grid<isochrone_padova> iso("Padova");
    for(double a=1.;a<9.;a+=0.02)
        std::cout<<a<<" "<<iso.max_mass(0.,a)<<" "<<iso.max_mass(-0.1,a)<<" "<<iso.max_mass(0.1,a)<<std::endl;
}
TEST(Gaia,Padova){
    isochrone_grid<isochrone_padova> iso("Padova", 1, 0.5);
    int NA=62; int NF=6;
    for(unsigned i=0;i<iso.iso(NF,NA)->N();++i)
        std::cout<<iso.iso(NF,NA)->mag(i,"G")<<" "
                 <<iso.iso(NF,NA)->mag(i,"K")<<std::endl;
}
TEST(MaxAge,Padova){
    isochrone_grid<isochrone_padova> iso("Padova", 1, 0.5);
    VecDoub metal_range = create_range(-3.,0.7,20);
    for(auto M:create_log_range(0.1,120.,30)){
        std::cout<<M<<" ";
        for(auto Z: metal_range)
            std::cout<<iso.max_age(Z,M)<<" ";
        std::cout<<std::endl;
    }
}
TEST(ES,Padova){
    isochrone_grid<isochrone_padova> iso("Padova",1,0.5);
    int NA=62; int NF=6;
    double age = iso.iso(NF,NA)->tau();
    double Z = iso.iso(NF,NA)->feh();
    std::cout<<Z<<" "<<age<<std::endl;
    double age1 = iso.iso(NF,NA+1)->tau();
    double Z1 = iso.iso(NF,NA+1)->feh();
    for(unsigned i=0;i<iso.iso(NF,NA)->N();++i){
        double mass = iso.iso(NF,NA)->initial_mass(i);
        double tmax = iso.max_age(Z,iso.iso(NF,NA)->initial_mass(i));
        std::cout<<mass<<" "<<age/tmax<<" "<<iso.iso(NF,NA)->logTeff(i)<<" "<<iso.iso(NF,NA)->logg(i)<<" ";
        mass = iso.iso(NF,NA+1)->initial_mass(i);
        tmax = iso.max_age(Z1,iso.iso(NF,NA+1)->initial_mass(i));
        std::cout<<mass<<" "<<age1/tmax<<" "<<iso.iso(NF,NA+1)->logTeff(i)<<" "<<iso.iso(NF,NA+1)->logg(i)<<" ";
        mass = iso.iso(NF,NA)->initial_mass(i);
        tmax = iso.max_age(Z1,iso.iso(NF,NA+1)->initial_mass(i));
        std::cout<<mass<<" "<<(age+0.25)/tmax<<" ";
        VecDoub ff = iso.interp_es(Z1,age+0.25,mass,{"V"});
        std::cout<<ff[0]<<" "<<ff[1]<<std::endl;
    }
}
TEST(ExtCurve,GExtCurve){
    schlafly2017_extinction_law EL;
    for(auto t: create_range(3.2,5.,10))   
        std::cout<<t<<" "<<EL.extinct_const("G",t)<<std::endl;
}
TEST(ExtCurve,ExtCurve){
    schlafly2017_extinction_law EL;
    std::vector<std::string> maglist = {"B","V","G","J","K","I","zP","Jv"};
    for(auto s: maglist)
        std::cout<<s<<" "<<EL.extinct_const(s)<<" "<<EL.extinct_const(s,3.6)<<" "<<EL.extinct_const(s,3.6,10.)<<std::endl;

}
TEST(ExtMap,ExtMapGrid){
    schlafly2017_extinction_law EL;
    sfd_extinction_map EM(&EL,conv::StandardSolarPAUL);
    for(double l=0.;l<2.*PI;l+=2.*PI/10.)
        for(double b=-PI/2.;b<PI/2.;b+=2.*PI/20.)
            for(double s = 0.01;s<2.;s+=0.1)
                std::cout<<l<<" "<<b<<" "<<s<<" "<<EM.A_V(l,b,s)<<std::endl;
}

TEST(ExtMapCombo,ExtMapGrid){
    schlafly2017_extinction_law EL;
    combo_extinction_map EM(&EL);
    for(double l=0.;l<2.*PI;l+=2.*PI/500.)
        for(double b=-PI/2.;b<PI/2.;b+=2.*PI/1000.)
            for(double s = 0.01;s<2.;s+=0.1)
                std::cout<<l<<" "<<b<<" "<<s<<" "<<EM.A_V(l,b,s)<<std::endl;
}

TEST(Params,Params){
    std::string p = parameters()["dir"]["extinction_coeffs"];
}

TEST(Prior, Prior){
    new_prior_2018 GP(conv::StandardSolarPAUL);
}

}
//=============================================================================

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
//=============================================================================
