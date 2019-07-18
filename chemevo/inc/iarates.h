#ifndef IARATES_H
#define IARATES_H
//=============================================================================
#include <map>
#include <string>
#include <cstddef>
//=============================================================================
#include "utils.h"
#include "params.h"
#include "imf.h"
#include "sfr.h"
#include "ages.h"
#include "radmig.h"
//=============================================================================
/**
 * @brief Supernova Type Ia Rate base class
 * @details Implements the base properties of the Supernova Type Ia rates. The
 * rates require an initial mass function, a star formation rate and a set of
 * stellar lifetimes. Also, to consider the migration of SN Type Ia progenitors
 * in the disc we also need a radial migration instance.
 */
class TypeIaRate{
private:
    std::shared_ptr<InitialMassFunction> imf;
    std::shared_ptr<StarFormationRate> sfr;
protected:
    std::shared_ptr<StellarLifetime> stellar_lifetime;
    std::shared_ptr<RadialMigration> radial_migration;
public:
    TypeIaRate(std::shared_ptr<InitialMassFunction> imf,
               std::shared_ptr<StarFormationRate> sfr,
               std::shared_ptr<StellarLifetime> SL,
               std::shared_ptr<RadialMigration> RM)
        :imf(imf),sfr(sfr),stellar_lifetime(SL),radial_migration(RM){}
    double SFR(double R, double t){return (*sfr)(R,t);}
    double IMF(double m){return (*imf)(m);}
    double RMKernel(double R, double Rp, double t){
        return (*radial_migration)(R,Rp,t);}
    virtual double operator()(double R, double t)=0;
};
//=============================================================================
/**
 * @brief Delay Time Distribution Supernova Type Ia Rate
 * @details Implements the Type Ia rate given by Greggio (2005) in terms of
 * a delay time distribution (DTD) that describes the probability of a Type Ia
 * dying at age tau. The rate is given by
 *
 * rate(R,t) = k_alpha \int_{tau_i}^{min(t,\tau_x)} A_B sfr(t-\tau)
 *                                                      DTD(\tau) d\tau
 *
 * where A_B is the fraction of all binary systems that make up the IMF
 * that give rise to Type Ia SN, \tau_i is the minimum delay time (set to
 * the lifetime of an 8 M_solar star) and \tau_x is the maximum delay time
 * (set to the lifetime of a 0.8 M_solar star).
 */
class TypeIaRate_DTD: public TypeIaRate{
protected:
    double ABinary = 0.0025;
    double tau_min = 0.04;
    double tau_max = 14.;
    double Norm;
public:
    TypeIaRate_DTD(ModelParameters M,
                   std::shared_ptr<InitialMassFunction> imf,
                   std::shared_ptr<StarFormationRate>sfr,
                   std::shared_ptr<RadialMigration> RM,
                   std::shared_ptr<StellarLifetime> SL=nullptr)
        :TypeIaRate(imf,sfr,SL,RM){
            ABinary=M.parameters["typeIa"]["BinaryFraction"];
            double MB_min=M.parameters["typeIa"]["MinimumIaBinaryMass"];
            double MB_max=M.parameters["typeIa"]["MaximumIaBinaryMass"];
            tau_min=(*SL)(MB_max,0.02);
            tau_max=(*SL)(MB_min,0.02);
            ABinary*=imf->norm();
        }
    /**
     * @brief delay time distribution
     * @details probability of a supernova type Ia dying with age \tau. Normalized such that \int_{tau_i}^{min(t,\tau_x)} DTD(\tau)d\tau = 1
     *
     * @param t age
     * @return probability of supernova dying at age tau.
     */
    virtual double DTD(double t)=0;
    /**
     * @brief Type Ia supernova rate
     *
     * @param R radius
     * @param t time
     *
     * @return number of Type Ia supernovae at radius R and time t
     */
    double operator()(double R, double t);
};
//=============================================================================
/**
 * @brief Matteucci et al. (2006) delay time distribution
 */
class MVP06_TypeIaRate:public TypeIaRate_DTD{
public:
    MVP06_TypeIaRate(ModelParameters M,
                     std::shared_ptr<InitialMassFunction> imf,
                     std::shared_ptr<StarFormationRate> sfr,
                     std::shared_ptr<RadialMigration> RM,
                     std::shared_ptr<StellarLifetime> SL=nullptr);
    double DTD(double t);
};
//=============================================================================
/**
 * @brief Supernova Type Ia rate in Matteucci & Greggio 1986 formulation.
 * @details Supernova Type Ia rate using the formulation from Matteucci &
 * Greggio 1986 where we explicitly consider the distribution of masses in
 * binary stars and count those that will form Type Ia supernovae.
 * The rate is given by
 *
 * Rate(R,t)=A_B\int_Ml^Mu dm imf(m)\int_{\mu_min}^0.5 f(mu) sfr(t-tau(m2)) dmu
 *
 * where A_B is the fraction of binaries in the mass range to form Type Ia SN
 * (between Mu=3 and Ml=16) that will form Type Ia SN and f(mu) is the
 * distribution of reduced masses given by
 *
 * f(mu) = 2^(1+gamma)(1+gamma)mu^(gamma).
 *
 * The more massive first mass m1 has evolved to a white dwarf and m2 is less
 * massive so stars overflowing its Roche lobe later (at age=tau(m2)) when the
 *  SN goes off. mu_min is the minimum mass fraction and is equal to the ratio
 *  of the mass of a star born at the earliest time which is dying now to the
 *  total mass of the binary. In some cases, this will mean that the mass of
 *  the primary is great enough that it has exploded as a Type II SN. In this
 *  case mu_min = (MB-min_mass_SNII)/MB.
 */
class TypeIaRate_BinaryMass: public TypeIaRate{
protected:
    double ABinary = 0.05;
    double MB_min = 3.;
    double MB_max = 16.;
    const double gamma = 2.;
    double tau_min;
public:
    TypeIaRate_BinaryMass(ModelParameters M,
                          std::shared_ptr<InitialMassFunction> imf,
                          std::shared_ptr<StarFormationRate> sfr,
                          std::shared_ptr<RadialMigration> RM,
                          std::shared_ptr<StellarLifetime> SL)
        :TypeIaRate(imf,sfr,SL,RM){
        ABinary=M.parameters["typeIa"]["BinaryFraction"];
        MB_min=M.parameters["typeIa"]["MinimumIaBinaryMass"];
        MB_max=M.parameters["typeIa"]["MaximumIaBinaryMass"];
        tau_min=(*SL)(MB_max,0.02);
    }
    double lifetime(double m){return (*stellar_lifetime)(m,0.02);}
    double min_mass_fraction(double t, double m);
    double f(double mu);
    double operator()(double R, double t);
};
//=============================================================================
extern unique_map<TypeIaRate,ModelParameters,
            std::shared_ptr<InitialMassFunction>,
            std::shared_ptr<StarFormationRate>,
            std::shared_ptr<RadialMigration>,
            std::shared_ptr<StellarLifetime>> tia_types;
//=============================================================================
/**
 * @brief Simple structure for computing the rate of type Ia SN in the
 * explicit binary formulation without radial migration.
 */
struct TypeIaRate_BinaryMass_st{
    TypeIaRate_BinaryMass *tia;
    double R;
    double t;
    double m;
};
/**
 * @brief Simple structure for computing the rate of type Ia SN in the
 * explicit binary formulation including radial migration (also pass birth
 * radius).
 */
struct TypeIaRate_BinaryMass_RM_st{
    TypeIaRate_BinaryMass *tia;
    double R;
    double t;
    double m;
    double Rp;
};
/**
 * @brief Simple structure for computing the rate of type Ia SN in the
 * explicit binary formulation including radial migration.
 */
struct TypeIaRate_BinaryMass_st_2D:TypeIaRate_BinaryMass_st{
    VecDoub x2min, x2max;
    TypeIaRate_BinaryMass_st_2D(TypeIaRate_BinaryMass *tia, double R, double t, double m, VecDoub x2min, VecDoub x2max):TypeIaRate_BinaryMass_st({tia,R,t,m}),x2min(x2min),x2max(x2max){}
};
//=============================================================================
/**
 * @brief Simple structure for computing the rate of type Ia SN in the
 * delay time distribution formulation without radial migration.
 */
struct TypeIaRate_DTD_st{
    TypeIaRate_DTD *tia;
    double R;
    double t;
};
/**
 * @brief Simple structure for computing the rate of type Ia SN in the
 * delay time distribution formulation including radial migration.
 */
struct TypeIaRate_DTD_st_2D:TypeIaRate_DTD_st{
    VecDoub x2min, x2max;
    TypeIaRate_DTD_st_2D(TypeIaRate_DTD *tia, double R, double t, VecDoub x2min, VecDoub x2max):TypeIaRate_DTD_st({tia,R,t}),x2min(x2min),x2max(x2max){}
};
//=============================================================================
#endif
//=============================================================================
