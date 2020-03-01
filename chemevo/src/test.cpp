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

TEST(StarFormationRate,SB15){
    SFR_SB15 sfr_sb15;
    double rt = sfr_sb15.gas_consumed_per_unit_radius(8.);
    double t = sfr_sb15.gas_consumed();
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
    ModelParameters M("example_params.json");
    M.print();
    M.pretty_print();
}

//=============================================================================

TEST(StellarAges,Portinari){

    ModelParameters MM("example_params.json");
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
    ModelParameters MM("example_params.json");
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
    ModelParameters MM("example_params.json");
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
    ModelParameters MM("example_params.json");
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
    ModelParameters M("example_params.json");
    NoneRadialMigration rm(M);
    MVP06_TypeIaRate tia(M,std::make_shared<ChabrierIMF>(cimf),
                              std::make_shared<SFR_SB15>(sfr_sb15),
                              std::make_shared<NoneRadialMigration>(rm));
    double t = tia(8.,5.);
}
TEST(TypeIaRates,BinaryMass){
    SFR_SB15 sfr_sb15;
    ChabrierIMF cimf;
    ModelParameters M("example_params.json");
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
    M.parameters["yields"]["AGB"]="Karakas";
    M.parameters["yields"]["typeIa"]="Maeda";
    M.parameters["yields"]["typeII"]="Kobayashi";
    YieldsSet Y(M);
    AGBYieldsKarakas agb(M);TypeIIKobayashi2006 typeII(M);TypeIaYields typeIa(M);
    EXPECT_NEAR(Y.mass("B",1.9,0.02),2.55e-25,1e-27);
    EXPECT_NEAR(Y.mass("Li",3.5,0.02),1.6e-10,1e-15);
    EXPECT_NEAR(Y.mass_remnant(1.,0.02),0.564,1e-10);
    EXPECT_NEAR(Y.mass_ejected(1.,0.02),1.-0.564,1e-10);
    EXPECT_NEAR(Y.mass("Ge",40.,0.02),0.002412,1e-14);
    EXPECT_NEAR(Y.mass("H",2.,0.01),agb.mass("H",2.,0.01),test_err);
    EXPECT_NEAR(Y.mass("O",20.,0.01),typeII.mass("O",20.,0.01),test_err);
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
    ModelParameters M("params/example_params.json");
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
    ModelParameters M("example_params.json");
    DoubleInfallInflow I(M);
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
    PezzulliInflowRadialFlow_rSFR<RadialFlow> P(M,sfr_sb15(8.3,13.7));
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
    ModelParameters M("example_params.json");
    std::shared_ptr<StarFormationRate> sfr = sfr_types["SB15"](M);
    std::shared_ptr<InitialMassFunction> imf = imf_types["Chabrier"](M);
    std::shared_ptr<StellarLifetime> sl = life_types["Portinari1998"](M);
    std::shared_ptr<RadialMigration> rm = rm_types["None"](M);
    std::unique_ptr<TypeIaRate> tia = tia_types["Matteucci2006"](M,imf,sfr,rm,sl);
}
TEST(Model, Model){
    ModelParameters P("example_params.json");
    Model M(P);
    P.parameters["fundamentals"]["IMF"]="WTF";
    try{Model M2(P);}
    catch(std::exception const & err){
        std::cerr<<err.what()<<std::endl;
        LOG(INFO)<<err.what()<<std::endl;
    }
}
TEST(Parameters,Parameters){
    ModelParameters P("example_params.json");
    try {double t = P.parameters["fundamentals"]["IM"];}
    catch(std::exception const & err){
        std::cerr<<err.what()<<std::endl;
        LOG(INFO)<<err.what()<<std::endl;
    }
}
TEST(SolarAbundance,SolarAbundance){
    ModelParameters P("example_params.json");
    AsplundSolarAbundances S(P);
    EXPECT_NEAR(S.Z(),0.0134,1e-4);
    AndersSolarAbundances S2(P);
    EXPECT_NEAR(S2.Z(),0.0189,1e-4);
}
//=============================================================================
TEST(RadialMigration, Rates){
    ModelParameters M("example_params.json");
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
    s2 = model4.DeathRate(8.3,12.);
    f2 = model4.EnrichmentRate("Fe",8.3,12.);
    f2 = model4.EnrichmentRate("H",8.3,12.);
    f2 = model4.GasReturnRate(0.3,.3);
    std::cout<<f2<<std::endl;
    a2 = model4.GasReturnRate(8.3,12.);
    sn2 = model4.SNIaRate(8.3,12.);
}
}
//=============================================================================

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  LOG(INFO)<<"Running tests";
  return RUN_ALL_TESTS();
}
//=============================================================================
