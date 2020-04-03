//=============================================================================
// Test programs
//=============================================================================
// #include "Python.h"
#include "utils_ch.h"
#include "utils.h"
#include "sfr.h"
#include "imf.h"
#include "ages.h"
#include "iarates.h"
#include "params.h"
#include "yields.h"
#include "grid.h"
#include "in_out.h"
#include "solar.h"
#include "model.h"
//=============================================================================
#include "gtest/gtest.h"
#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

const double test_err = 1e-7;

//=============================================================================

namespace {
//=============================================================================
TEST(GasDump,GasDump){
    ModelParameters M("for_kai/kai_test_double_infall.json");
    std::shared_ptr<SolarAbundances> S = solar_types["Asplund"](M);
    std::shared_ptr<StarFormationRate> SFR = sfr_types[M.parameters["fundamentals"]["SFR"]](M);
    GasDumpSimple gds(M,S);
    for(auto r: create_range(0.1,13.,50))
        std::cout<<r<< " "<< (*SFR)(8.3,r)<<std::endl;
}
TEST(StarFormationRate,DoubleInfall){
    ModelParameters M("params/example_params.json");
    SFR_DoubleInfall di(M);
    std::cout<<di(0.2,0.)<<std::endl;
    std::cout<<di(8.2,10.)<<std::endl;
    std::cout<<di(8.2,12.)<<std::endl;
}
TEST(StarFormationRate,SB15){
    SFR_SB15 sfr_sb15;
    double rt = sfr_sb15.gas_consumed_per_unit_radius(8.);
    double t = sfr_sb15.gas_consumed();
}
TEST(StarFormationRate,Neige2020){
    ModelParameters M("params/example_params.json");
    SFR_Neige2020 neige20(M);
    std::cout<<neige20(5.2,10.)<<std::endl;
    std::cout<<neige20(8.2,10.)<<std::endl;
    std::cout<<neige20(8.2,12.)<<std::endl;
    double rt = neige20.gas_consumed_per_unit_radius(8.);
    double t = neige20.gas_consumed();
}
//=============================================================================

TEST(InitialMassFunction,Salpeter){
    SalpeterIMF kimf;
    double rf = kimf(1.);
    double t = kimf.mean_mass();
    EXPECT_NEAR(t,1.,1e-4);
}
TEST(InitialMassFunction,Tinsley){
    TinsleyIMF kimf;
    double rf = kimf(1.);
    double t = kimf.mean_mass();
    EXPECT_NEAR(t,1.,1e-4);
}
TEST(InitialMassFunction,Scalo){
    ScaloIMF kimf;
    double rf = kimf(1.);
    double t = kimf.mean_mass();
    EXPECT_NEAR(t,1.,1e-4);
}
TEST(InitialMassFunction,Kroupa){
    KroupaIMF kimf;
    double rf = kimf(1.);
    double t = kimf.mean_mass();
    EXPECT_NEAR(t,1.,1e-4);
}
TEST(InitialMassFunction,Chabrier){
    ChabrierIMF kimf;
    double rf = kimf(1.);
    double t = kimf.mean_mass();
    EXPECT_NEAR(t,1.,1e-4);
}
//=============================================================================

TEST(Galaxy,grid){
    Grid grid(10,40,0.1,10.,12.);
    SFR_SB15 sfr_sb15;
    VecDoub rr = grid.grid_radial();
    VecDoub tt = grid.grid_time();
    for(unsigned t=0;t<tt.size();++t){
	    for(unsigned r=0;r<rr.size();++r)
	    	grid.set(sfr_sb15(rr[r],tt[t]),r,t);
	}
}
//=============================================================================
TEST(Parameters,params){
    ModelParameters M("params/example_params.json");
    M.print();
    M.pretty_print();
}
//=============================================================================
TEST(TypeIaYields, Maeda){
    ModelParameters M("params/example_params.json");
    Maeda_TypeIaYields I(M);
    EXPECT_NEAR(I.mass("C"),0.0499,1e-7);
}
TEST(TypeIaYields, Iwamoto){
    ModelParameters M("params/example_params.json");
    Iwamoto_TypeIaYields I(M);
    EXPECT_NEAR(I.mass("C"),0.0483014,1e-7);
}
TEST(TypeIaYields, Seitenzahl){
    ModelParameters M("params/example_params.json");
    Seitenzahl_TypeIaYields I(M);
    EXPECT_NEAR(I.mass("C"),3.04e-3,1e-5);
}
TEST(TypeIaYields, Thielemann){
    ModelParameters M("params/example_params.json");
    Thielemann_TypeIaYields I(M);
    EXPECT_NEAR(I.mass("C"),0.0504011,1e-7);
}
TEST(TypeIIYields, Nugrid){
    ModelParameters M("params/example_params.json");
    M.parameters["fundamentals"]["solar"]="GrevesseNoel";
    NugridYields I(M);
    EXPECT_NEAR(I.mass("Ce",20.,0.02),9.710E-08,1e-12);
    VecDoub mm = create_range(1.,50.,50);
    for(auto i: mm)
        std::cout<<i<<" "<<I.mass("Fe",i,0.001)<<std::endl;
}
TEST(TypeIIYields, Kobayashi2006){
    ModelParameters M("params/example_params.json");
    M.parameters["fundamentals"]["solar"]="Anders";
    TypeIIKobayashi2006 I(M);
    EXPECT_NEAR(I.mass("C",15.,0.02),0.06625,1e-5);
    VecDoub mm = create_range(8.,50.,50);
    for(auto i: mm)
        std::cout<<i<<" "<<I.mass("Fe",i,0.02)<<std::endl;
}
TEST(TypeIIYields, CompareYields){
    ModelParameters M("params/example_params.json");
    M.parameters["fundamentals"]["solar"]="Asplund";
    TypeIIKobayashi2006 K(M);
    TypeIIChieffiLimongi2004 CL04(M);
    TypeIIChieffiLimongi2018 CL18(M);
    NugridYields N(M);
    VecDoub mm = create_range(8.,50.,50);
    std::string el = "Fe";
    double Z = pow(10,-1)*0.0198;
    for(auto i: mm)
        std::cout<<i<<" "
        <<K.mass(el,i,Z)<<" "
        <<CL04.mass(el,i,Z)<<" "
        <<CL18.mass(el,i,Z)<<" "
        <<N.mass(el,i,Z)<<
        std::endl;
}
TEST(AGBYields, CompareYields){
    ModelParameters M("params/example_params.json");
    M.parameters["fundamentals"]["solar"]="Asplund";
    AGBYieldsKarakas K(M);
    AGBYieldsVentura V(M);
    VecDoub mm = create_range(1.,6.,50);
    std::string el = "He";
    double Z = 0.02;
    for(auto i: mm)
        std::cout<<i<<" "
        <<K.mass(el,i,Z)<<" "
        <<V.mass(el,i,Z)<<
        std::endl;
}
TEST(TypeIIYields, ChieffiLimongi2004){
    ModelParameters M("params/example_params.json");
    M.parameters["fundamentals"]["solar"]="Anders";
    TypeIIChieffiLimongi2004 I(M);
    EXPECT_NEAR(I.mass("C",15.,0.02),0.20631,1e-5);
    EXPECT_NEAR(I.yield("C",15.,0.02),0.0096337751473700058,1e-5);
}
TEST(TypeIIYields, ChieffiLimongi2018){
    ModelParameters M("params/example_params.json");
    M.parameters["fundamentals"]["solar"]="Asplund";
    TypeIIChieffiLimongi2018 I(M);
    EXPECT_NEAR(I.mass("C",25.,1.345e-2),0.620690828,1e-5);
}
TEST(AGBYields, KarakasYields){
    ModelParameters M("params/example_params.json");
    M.parameters["fundamentals"]["solar"]="Anders";
    AGBYieldsKarakas I(M);
    EXPECT_NEAR(I.mass("C",1.,0.02),1.1907e-3,1e-7);
    EXPECT_NEAR(I.yield("C",1.,0.02),-0.00082472786359013965,1e-7);
}
TEST(AGBYields, VenturaYields){
    ModelParameters M("params/example_params.json");
    M.parameters["fundamentals"]["solar"]="Asplund";//GrevesseSauval";
    AGBYieldsVentura I(M);
    EXPECT_NEAR(I.mass("H",2.5,0.02),7.32e-7,1e-9);
}
//=============================================================================

TEST(StellarAges,Portinari){

    ModelParameters MM("params/example_params.json");
    Portinari1998 lifetime(MM);
    EXPECT_NEAR(lifetime(1.1,0.004),4.93,1e-3);
    EXPECT_NEAR(lifetime(40.,0.02),5.12/1000.,1e-8);
    double tmp = lifetime(150.,0.001);

    EXPECT_NEAR(lifetime.mass_star_dying_now(lifetime(40.,0.02),0.02),40.,5e-3);
    EXPECT_NEAR(lifetime.mass_star_dying_now(lifetime(40.,0.012),0.012),40.,0.1);
    EXPECT_NEAR(lifetime.mass_star_dying_now(lifetime(4.,0.002),0.002),4.,0.05);
    double Mup=lifetime.mass_star_dying_now(4.1,0.002);
    double Mdo=lifetime.mass_star_dying_now(3.9,0.002);
    double M=lifetime.mass_star_dying_now(4.,0.002);
    EXPECT_NEAR(lifetime.dlogMdlogt(4.,0.002)*M/4.,(Mup-Mdo)/0.2,1e-4);
    Mup=lifetime.mass_star_dying_now(0.0201,0.002);
    Mdo=lifetime.mass_star_dying_now(0.0199,0.002);
    M=lifetime.mass_star_dying_now(0.02,0.002);
    EXPECT_NEAR(lifetime.dlogMdlogt(0.02,0.002)*M/0.02,(Mup-Mdo)/0.0002,2e-2);

}
TEST(StellarAges,PadovaniMatteucci){
    ModelParameters MM("params/example_params.json");
    PadovaniMatteucci1993 lifetime(MM);

    EXPECT_NEAR(lifetime.mass_star_dying_now(lifetime(40.,0.02),0.02),40.,1e-2);
    EXPECT_NEAR(lifetime.mass_star_dying_now(lifetime(4.,0.002),0.002),4.,1e-2);
    double Mup=lifetime.mass_star_dying_now(4.1,0.002);
    double Mdo=lifetime.mass_star_dying_now(3.9,0.002);
    double M=lifetime.mass_star_dying_now(4.,0.002);
    EXPECT_NEAR(lifetime.dlogMdlogt(4.,0.002)*M/4.,(Mup-Mdo)/0.2,1e-4);
    Mup=lifetime.mass_star_dying_now(0.0201,0.002);
    Mdo=lifetime.mass_star_dying_now(0.0199,0.002);
    M=lifetime.mass_star_dying_now(0.02,0.002);
    EXPECT_NEAR(lifetime.dlogMdlogt(0.02,0.002)*M/0.02,(Mup-Mdo)/0.0002,1e-2);

}
TEST(StellarAges,Kodama1997){
    ModelParameters MM("params/example_params.json");
    Kodama1997 lifetime(MM);

    EXPECT_NEAR(lifetime.mass_star_dying_now(lifetime(40.,0.02),0.02),40.,1e-2);
    EXPECT_NEAR(lifetime.mass_star_dying_now(lifetime(4.,0.002),0.002),4.,1e-2);
    double Mup=lifetime.mass_star_dying_now(4.1,0.002);
    double Mdo=lifetime.mass_star_dying_now(3.9,0.002);
    double M=lifetime.mass_star_dying_now(4.,0.002);
    EXPECT_NEAR(lifetime.dlogMdlogt(4.,0.002)*M/4.,(Mup-Mdo)/0.2,1e-4);
    Mup=lifetime.mass_star_dying_now(0.0201,0.002);
    Mdo=lifetime.mass_star_dying_now(0.0199,0.002);
    M=lifetime.mass_star_dying_now(0.02,0.002);
    EXPECT_NEAR(lifetime.dlogMdlogt(0.02,0.002)*M/0.02,(Mup-Mdo)/0.0002,1e-2);
}
TEST(StellarAges,MaederMeynet1989){
    ModelParameters MM("params/example_params.json");
    MaederMeynet1989 lifetime(MM);

    EXPECT_NEAR(lifetime.mass_star_dying_now(lifetime(40.,0.02),0.02),40.,1e-2);
    EXPECT_NEAR(lifetime.mass_star_dying_now(lifetime(4.,0.002),0.002),4.,1e-2);
    double Mup=lifetime.mass_star_dying_now(4.1,0.002);
    double Mdo=lifetime.mass_star_dying_now(3.9,0.002);
    double M=lifetime.mass_star_dying_now(4.,0.002);
    EXPECT_NEAR(lifetime.dlogMdlogt(4.,0.002)*M/4.,(Mup-Mdo)/0.2,1e-4);
    Mup=lifetime.mass_star_dying_now(0.0201,0.002);
    Mdo=lifetime.mass_star_dying_now(0.0199,0.002);
    M=lifetime.mass_star_dying_now(0.02,0.002);
    EXPECT_NEAR(lifetime.dlogMdlogt(0.02,0.002)*M/0.02,(Mup-Mdo)/0.0002,1e-2);
}
//=============================================================================

TEST(TypeIaRates,MVP06){
    SFR_SB15 sfr_sb15;
    ChabrierIMF cimf;
    ModelParameters M("params/example_params.json");
    NoneRadialMigration rm(M);
    MVP06_TypeIaRate tia(M,std::make_shared<ChabrierIMF>(cimf),
                              std::make_shared<SFR_SB15>(sfr_sb15),
                              std::make_shared<NoneRadialMigration>(rm));
    double t = tia(8.,5.);
}
TEST(TypeIaRates,BinaryMass){
    SFR_SB15 sfr_sb15;
    ChabrierIMF cimf;
    ModelParameters M("params/example_params.json");
    NoneRadialMigration rm(M);
    Portinari1998 lifetime(M);
    TypeIaRate_BinaryMass tia(M,std::make_shared<ChabrierIMF>(cimf),
                              std::make_shared<SFR_SB15>(sfr_sb15),
                              std::make_shared<NoneRadialMigration>(rm),
                              std::make_shared<Portinari1998>(lifetime));
    double t = tia(8.,5.);
}
//=============================================================================

TEST(Yields,YieldsSet){
    ModelParameters M("params/example_params.json");
    M.parameters["yields"]["AGB"]="Nugrid";
    M.parameters["yields"]["typeIa"]="Maeda";
    M.parameters["yields"]["typeII"]="Kobayashi";
    YieldsSet Y(M);exit(1);
    AGBYieldsKarakas agb(M);
    TypeIIKobayashi2006 typeII(M);
    Maeda_TypeIaYields typeIa(M);
    EXPECT_NEAR(Y.mass("B",1.9,0.02),2.55e-25,1e-27);
    EXPECT_NEAR(Y.mass("Li",3.5,0.02),1.6e-10,1e-15);
    EXPECT_NEAR(Y.mass_remnant(1.,0.02),0.564,1e-10);
    EXPECT_NEAR(Y.mass_ejected(1.,0.02),1.-0.564,1e-10);
    EXPECT_NEAR(Y.mass("Ge",40.,0.02),0.002412,1e-14);
    EXPECT_NEAR(Y.mass("H",2.,0.01),agb.mass("H",2.,0.01),test_err);
    EXPECT_NEAR(Y.mass("O",20.,0.01),typeII.mass("O",20.,0.01),test_err);
    std::cout<<Y.typeIa_ejectedmass("Ni")<<std::endl;
    EXPECT_NEAR(Y.typeIa_ejectedmass("Ni"),typeIa.mass("Ni"),test_err);
    for(double mm=1.;mm<8.;mm+=0.4)
        std::cout<<mm<<" "<<Y.mass_ejected(mm,0.001)<<" "<<Y.mass("H",mm,0.001)<<std::endl;
}
TEST(Yields,ChieffiLimongi){
    ModelParameters M("params/example_params.json");
    M.parameters["yields"]["typeII"]="ChieffiLimongi";
    M.parameters["yields"]["AGB"]="Karakas";
    M.parameters["yields"]["typeIa"]="Maeda";
    YieldsSet Y(M);
    EXPECT_NEAR(Y.mass("F",15.,0.02),4.19E-06,1e-13);
}

TEST(Yields,AGBYields){
    ModelParameters M("params/params/example_params.json");
    M.parameters["yields"]["AGB"]="Karakas";
    AGBYieldsKarakas agb(M);
    for(double mm=0.4;mm<8.;mm+=0.2)
        std::cerr<<mm<<" "
                 <<agb.mass_ejected(mm,0.001)<<" "
                 <<agb.mass_ejected(mm,0.01)<<" "
                 <<agb.mass_ejected(mm,0.1)<<" "
                 <<agb.mass("He",mm,0.001)<<" "
                 <<agb.mass("He",mm,0.01)<<" "
                 <<agb.mass("He",mm,0.1)<<std::endl;
}
//=============================================================================

TEST(Flows, flows){
    ModelParameters M("params/example_params.json");
    std::shared_ptr<AndersSolarAbundances> solar;
    DoubleInfallInflow I(M, solar);
    SimpleGalacticFountain G(M);
    LinearRadialFlow L(M);
}

// TEST(Flows, Pezzulliflows){
//     ModelParameters M("example_params_pezzulli.json");
//     std::shared_ptr<StarFormationRate> sfr = sfr_types["ExpDecay"](M);
//     PezzulliInflowRadialFlow P(M,sfr);
//     EXPECT_NEAR(P.mu(14.,13.7),-83.34095767,1e-3);
//     EXPECT_NEAR(P.mu(2.,13.7),-124.134584719,1e-3);
//     EXPECT_NEAR(P.acc_rate(4.,13.7),8.0551719,1e-4);
//     EXPECT_NEAR(P.acc_rate(10.,13.7),0.94847,1e-4);
//     EXPECT_NEAR(P.flow_rate(2.,4.),-0.192824798,1e-4);
//     EXPECT_NEAR(P.flow_rate(8.,10.),-0.40291924,1e-4);
// }
TEST(Flows, PezzulliReducedSFRFlows){
    ModelParameters M("example_params_pezzulli.json");
    Grid grid(M);
    SFR_ExpDecay sfr_sb15(M);
    VecDoub rr = grid.grid_radial();
    VecDoub tt = grid.grid_time();
    for(unsigned t=0;t<tt.size();++t){
        for(unsigned r=0;r<rr.size();++r)
            grid.set(sfr_sb15(rr[r],tt[t]),r,t);
    }
    std::shared_ptr<AndersSolarAbundances> solar;
    PezzulliInflowRadialFlow_rSFR<RadialFlow> P(M,solar,sfr_sb15(8.3,13.7));
    std::shared_ptr<StarFormationRate> sfr = sfr_types["ExpDecay"](M);
    // PezzulliInflowRadialFlow P2(M,sfr);

    // EXPECT_NEAR(P.sigmagas(8.,5.,&grid),P2.sigmagas(8.,5.),1e-2);
    // EXPECT_NEAR(P.sigmaeffdot(8.,5.,&grid),P2.sigmaeffdot(8.,5.),1e-2);

    EXPECT_NEAR(P.mu(14.,13.7,&grid),-83.34095767,1.);
    EXPECT_NEAR(P.mu(2.,13.7,&grid),-124.134584719,1.);
    EXPECT_NEAR(P.acc_rate(4.,13.7,&grid),8.0551719,3e-2);
    EXPECT_NEAR(P.acc_rate(10.,13.7,&grid),0.94847,3e-2);
    EXPECT_NEAR(P.flow_rate(2.,4.,&grid),-0.192824798,3e-2);
    EXPECT_NEAR(P.flow_rate(8.,10.,&grid),-0.40291924,3e-2);
    // EXPECT_NEAR(P.acc_rate(20.,13.7,&grid),P2.acc_rate(20.,13.7),3e-3);
    // EXPECT_NEAR(P.flow_rate(20.,13.7,&grid),P2.flow_rate(20.,13.7),3e-2);

}
//=============================================================================

TEST(Maps, maps){
    ModelParameters M("params/example_params.json");
    std::shared_ptr<StarFormationRate> sfr = sfr_types["SB15"](M);
    std::shared_ptr<InitialMassFunction> imf = imf_types["Chabrier"](M);
    std::shared_ptr<StellarLifetime> sl = life_types["Portinari1998"](M);
    std::shared_ptr<RadialMigration> rm = rm_types["None"](M);
    std::unique_ptr<TypeIaRate> tia = tia_types["Matteucci2006"](M,imf,sfr,rm,sl);
}
TEST(Model, Model){
    ModelParameters P("params/example_params.json");
    Model M(P);
    P.parameters["fundamentals"]["IMF"]="WTF";
    try{Model M2(P);}
    catch(std::exception const & err){
        std::cerr<<err.what()<<std::endl;
        LOG(INFO)<<err.what()<<std::endl;
    }
}
TEST(Parameters,Parameters){
    ModelParameters P("params/example_params.json");
    try {double t = P.parameters["fundamentals"]["IM"];}
    catch(std::exception const & err){
        std::cerr<<err.what()<<std::endl;
        LOG(INFO)<<err.what()<<std::endl;
    }
}
TEST(SolarAbundance,SolarAbundance){
    ModelParameters P("params/example_params.json");
    AsplundSolarAbundances S(P);
    EXPECT_NEAR(S.Z(),0.0134,1e-4);
    AndersSolarAbundances S2(P);
    EXPECT_NEAR(S2.Z(),0.0189,1e-4);
}


//=============================================================================
TEST(RadialMigration, BasicTest){
    ModelParameters M("params/example_params.json");
    M.parameters["migration"]["Form"]="GaussianDrift";
    M.parameters["migration"]["sigmaR"]=4.0;
    M.parameters["fundamentals"]["StarScaleLength"]=1.;
    M.parameters["fundamentals"]["GalaxyAge"]=1.;
    Grid gas_mass(50,3,0.1,20.,2.);
    Grid mass_fraction(50,3,0.1,20.,2.);
    for(auto i=0;i<50;++i){
        gas_mass.set(exp(-gas_mass.grid_radial()[i]),i,0);
        mass_fraction.set(0.01,i,0);
    }
    unsigned nRR = 10, nt=1;
    std::shared_ptr<RadialMigration> RM = rm_types[M.parameters["migration"]["Form"]](M);
    for(unsigned nRR=0;nRR<50;++nRR)
        std::cout<<gas_mass.grid_time()[nt]<<" "<<gas_mass.grid_radial()[nRR]<<" "
        <<RM->convolve(&gas_mass,nRR,nt)-gas_mass(nRR,nt-1)<<" "<<RM->convolve_massfrac(&gas_mass,&mass_fraction,nRR,nt)<<std::endl;
}
//=============================================================================
TEST(RadialMigration, Rates){
    ModelParameters M("params/example_params.json");
    M.parameters["yields"]["typeII"]="ChieffiLimongi2004";
    M.parameters["yields"]["AGB"]="Karakas";
    M.parameters["yields"]["typeIa"]="Maeda";
    M.parameters["migration"]["Form"]="None";
    Model model(M);
    auto s1 = model.DeathRate(8.3,12.);
    auto f1 = model.EnrichmentRate("Fe",8.3,12.);
    auto a1 = model.GasReturnRate(8.3,12.);
    auto sn1 = model.SNIaRate(8.3,12.);
    M.parameters["migration"]["Form"]="Gaussian";
    M.parameters["migration"]["sigmaR"]=0.05;
    Model model2(M);
    auto s2 = model2.DeathRate(8.3,12.);
    auto f2 = model2.EnrichmentRate("Fe",8.3,12.);
    auto a2 = model2.GasReturnRate(8.3,12.);
    auto sn2 = model2.SNIaRate(8.3,12.);
    EXPECT_NEAR(s1,s2,3e-3);
    EXPECT_NEAR(f1,f2,1e-4);
    EXPECT_NEAR(a1,a2,1e-3);
    EXPECT_NEAR(sn1,sn2,2e-5);
    M.parameters["migration"]["Form"]="GaussianDrift";
    M.parameters["migration"]["sigmaR"]=0.05;
    Model model3(M);
    s2 = model3.DeathRate(8.3,12.);
    f2 = model3.EnrichmentRate("Fe",8.3,12.);
    a2 = model3.GasReturnRate(8.3,12.);
    sn2 = model3.SNIaRate(8.3,12.);
    EXPECT_NEAR(s1,s2,3e-3);
    EXPECT_NEAR(f1,f2,1e-4);
    EXPECT_NEAR(a1,a2,1e-3);
    EXPECT_NEAR(sn1,sn2,2e-5);

    M.parameters["migration"]["Form"]="GaussianDrift";
    M.parameters["migration"]["sigmaR"]=5.0;
    Model model4(M);
    double R = 25.1;
    s2 = model4.DeathRate(R,12.);
    f2 = model4.EnrichmentRate("Fe",R,12.);
    f2 = model4.EnrichmentRate("H",R,12.);
    f2 = model4.GasReturnRate(R,.3);
    std::cout<<f2<<std::endl;
    a2 = model4.GasReturnRate(R,12.);
    sn2 = model4.SNIaRate(R,12.);
}
}
//=============================================================================

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  LOG(INFO)<<"Running tests";
  return RUN_ALL_TESTS();
}
//=============================================================================
