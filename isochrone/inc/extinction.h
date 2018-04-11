#ifndef EXT_H
#define EXT_H
//=============================================================================
#include "utils.h"
#include "utils_iso.h"
#include "coordtransforms.h"
//=============================================================================
#include "json.hpp"
#include "healpix_base.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "fitshandle.h"
//=============================================================================

class extinction_law{
private:
    //=========================================================================
    const double ABconst=1.324,AVconst=1.,AIconst=0.482,AJconst=0.282;
    const double AHconst=0.175,AKconst=0.112;
    const double AGconst=1.; // BAD
protected:
    std::map<std::string, double> coeff_list, avgradient;
public:
    extinction_law(void);
    double extinct_const(std::string band, double lTeff=4., double AV=0.);
    double av_gradient(std::string band);
    double extinct_const_colour(std::vector<std::string> colour,
                                double lTeff=4., double AV=0.);
    virtual double G_extinction(double lTeff, double AV=0.);
};

class schlafly2017_extinction_law: public extinction_law{
private:
    //========================================================================
    VecDoub G_lTeff, G_extinct; // For G band extinction
    const std::string extinction_coeff_file = "extinction_coeffs_2017.dat";
    const std::string G_extinction_coeff_file = "extinction_coeffs_G_2017.dat";
    const std::string avgradient_extinction_coeff_file = "avgradient_2017.dat";
public:
    //========================================================================
    schlafly2017_extinction_law(double RV=3.3);
    double G_extinction(double lTeff, double AV=0.);
};

//=============================================================================
/**
 * @brief Class to hold an extinction map
 */

class extinction_map{
private:
    extinction_law *EL;           // Extinction law
public:
    extinction_map(extinction_law *EL):EL(EL){}
    double extinct_const(std::string band, double lTeff=4., double AV=0.){
        return EL->extinct_const(band,lTeff,AV);
    }
    double extinct_const_colour(std::vector<std::string> colour,
                                double lTeff=4., double AV=0.){
        return EL->extinct_const_colour(colour, lTeff, AV);
    }
    //========================================================================
    // Getters
    virtual double A_V(double l, double b, double s) = 0;

    inline double A_B(double l, double b, double s){
        return EL->extinct_const("B")*A_V(l,b,s);
    }
    inline double A_I(double l, double b, double s){
        return EL->extinct_const("I")*A_V(l,b,s);
    }
    inline double A_J(double l, double b, double s){
        return EL->extinct_const("J")*A_V(l,b,s);
    }
    inline double A_H(double l, double b, double s){
        return EL->extinct_const("H")*A_V(l,b,s);
    }
    inline double A_K(double l, double b, double s){
        return EL->extinct_const("K")*A_V(l,b,s);
    }
    inline double A_JK(double l, double b, double s){
        return EL->extinct_const_colour({"J","K"})*A_V(l,b,s);
    }
};

class sfd_extinction_map: public extinction_map{
private:
    const int Nl = 360, Nb = 180; // Resolution of grids in Galactic coords
    const int Ns = 50;            // Resolution of grid in log distance
    VecDoub LL, BB, SS;           // Grids in Galactic coords and log distance
    std::vector<std::vector<VecDoub>> MAP; // Extinction as a function of l,b,s
    bool use_factor;              // Use adjustment of Schlegel
    // Advice -- use RV=2.742 (Schlafly & Finkbeiner 2011)
    // or use RV=3.1 with use_factor=True
    //========================================================================
    const int NSIDE = 512;        // Property of healpix map
    // path to healpix map
    const std::string filename = parameters()["dir"]["extinction"];
    // "/data/jls/extinction/lambda_sfd_ebv.fits";
    Healpix_Map<double> HP_Map;   // Healpix extinction map
    //========================================================================
    // Parameters of global gas distribution
    const double hR = 4.2, hz=0.088, gammafl=0.0054,Rw=8.4,gammaw=0.18;
    double Rfl;
    VecDoub StandardSolar;
    double RV; // Extinction constant (conversion from Schlegel to AV)
    //========================================================================
    /**
     * @brief compute extinction in A_V
     * @details A_{Vinf} \int_0^s \d s rho(s)/ \int_0^\inf \d s rho(s)
     *          where A_{Vinf} is from Schlegel, Finkbeiner & Davis (1998)
     *
     * @param l Galactic l
     * @param b Galactic b
     * @param s distance
     * @return extinction in V band
     */
    double evaluate(double l, double b, double s);
    /**
     * @brief factor by which Schlegel maps are modified
     *
     * @param EBV extinction in (B-V)
     */
    double factor(double EBV);
public:
    //========================================================================
    sfd_extinction_map(extinction_law *EL, VecDoub StS=conv::StandardSolarPAUL,
                   double RV=2.742, bool use_factor=false);
    //========================================================================
    /**
     * @brief global 3D density distribution of the gas
     *
     * @param l Galactic l
     * @param b Galactic b
     * @param s distance
     * @return gas density
     */
    double density(double l, double b, double s);

    double A_V(double l, double b, double s){
        return evaluate(l,b,s);
    }
};

class combo_extinction_map: public extinction_map{
private:
    const int NSIDE = 1024;                   // Healpix resolution
    const int npix  = 12*NSIDE*NSIDE;       // No. of healpix
    const int Nd    = 31;                   // Number of distance bins
    std::vector<Healpix_Map<double>> HP_Map_grid;
    VecDoub log10distance;                     // ln(distance) grid
    const std::string filename = parameters()["dir"]["extinction"];
    double evaluate(double l, double b, double s, bool interp=false);
public:
    combo_extinction_map(extinction_law *EL);
    double A_V(double l, double b, double s){
        return evaluate(l, b, s);
    }
};

#endif
