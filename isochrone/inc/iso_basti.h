#ifndef ISOBASTI_H
#define ISOBASTI_H
//=============================================================================
#include "iso_base.h"
//=============================================================================
// Data for BaSTI list of isochrones
namespace BaSTI{
    const VecDoub ZList = {0.00001,0.0001,0.0003,0.0006,0.001,0.002,0.004,0.008,0.01,0.0198,0.03,0.04};
    const VecDoub FeHList = {-3.27,-2.27,-1.79,-1.49,-1.27,-0.96,-0.66,-0.35,-0.25,0.06,  0.26,0.4};
    const VecDoub ageList = {0.05,0.10,0.5,1.0,2.0,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.};
}
namespace BaSTI_More{
    const VecDoub ZList = {0.00001,0.0001,0.0003,0.0006,0.001,0.002,0.004,0.008,0.01,0.0198,0.03,0.04};
    const VecDoub FeHList = {-3.27,-2.27,-1.79,-1.49,-1.27,-0.96,-0.66,-0.35,-0.25,0.06,  0.26,0.4};
    const VecDoub ageList = {0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12.};
}

//=============================================================================
/**
 * @brief Class to store a single BaSTI Johnson bands (BVIJHK) isochrone
 */
class isochrone_johnson: public isochrone{
    public:
       isochrone_johnson(void);

       double get_metallicity(std::vector<std::string> input_iso, std::string dir);
        void fill(std::vector<std::string> input_iso, std::string dir,
                  double Age=0.,double thin_mag=-1.);
};
//=============================================================================
#endif
//=============================================================================
