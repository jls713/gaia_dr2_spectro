#ifndef SF_IMPL_H
#define SF_IMPL_H

#include "utils_iso.h"
#include "sf.h"

extern std::string jason_dir(void);
extern std::string gus_dir(void);
extern std::string rave_dir(void);

extern std::map<std::string, Selection_function*(*)(void)> sf_types;
template<typename T> Selection_function * createSFInstance(void) {
    return new T();}

class RAVE_selection_function:public Binned_selection_function{
public:
    RAVE_selection_function(void):Binned_selection_function(69,0.,13.8,"I",32,rave_dir()+"RAVE_completeness.fits"){}
};

class TGAS_selection_function: public Binned_selection_function{
public:
    TGAS_selection_function(void):Binned_selection_function(260,2.,15.,"Vp",8,jason_dir()+"NSIDE8_map_Vp.fits"){}
};

class TGAS_HQ_selection_function: public Binned_selection_function{
public:
    TGAS_HQ_selection_function(void):Binned_selection_function(260,2.,15.,"Vp",8,jason_dir()+"NSIDE8_HQ_map_Vp.fits"){}
};

class Tycho2_selection_function8:public Binned_selection_function{
public:
    Tycho2_selection_function8(void):Binned_selection_function(30,4.,13.,"J",8,gus_dir()+"tgas_nside8_map.fits"){}
};

class Tycho2_selection_function32:public Binned_selection_function{
public:
    Tycho2_selection_function32(void):Binned_selection_function(30,4.,13.,"J",32,gus_dir()+"tgas_nside32_map.fits"){}
};

class Tycho2_selection_function:public Hybrid_Binned_selection_function{
public:
    Tycho2_selection_function(void):Hybrid_Binned_selection_function(8.,{30,30},{4.,4.},{13.,13.},"J",{8,32},{gus_dir()+"tgas_nside8_map.fits",gus_dir()+"tgas_nside32_map.fits"}){}
};

class Tycho2RAVE_selection_function8:public Binned_selection_function{
public:
    Tycho2RAVE_selection_function8(void):Binned_selection_function(30,2.56,12.18,"J",8,gus_dir()+"tycho_rave_nside8.fits"){}
};

class Tycho2RAVE_selection_function32:public Binned_selection_function{
public:
    Tycho2RAVE_selection_function32(void):Binned_selection_function(30,2.56,12.18,"J",32,gus_dir()+"tycho_rave_nside32.fits"){}
};

class Tycho2RAVE_selection_function:public Hybrid_Binned_selection_function{
public:
    Tycho2RAVE_selection_function(void):Hybrid_Binned_selection_function(8.,{30,30},{2.56,2.56},{12.18,12.18},"J",{8,32},{gus_dir()+"tycho_rave_nside8.fits",gus_dir()+"tycho_rave_nside32.fits"}){}
};

class Tycho2RAVE_selection_function8_I:public Binned_selection_function{
public:
    Tycho2RAVE_selection_function8_I(void):Binned_selection_function(30,6.,12.,"I",8,gus_dir()+"tycho_rave_nside8_Iband.fits"){}
};

class Tycho2RAVE_selection_function32_I:public Binned_selection_function{
public:
    Tycho2RAVE_selection_function32_I(void):Binned_selection_function(30,6.,12.,"I",32,gus_dir()+"tycho_rave_nside32_Iband.fits"){}
};

class Tycho2RAVE_selection_function_I:public Hybrid_Binned_selection_function{
public:
    Tycho2RAVE_selection_function_I(void):Hybrid_Binned_selection_function(8.,{30,30},{6.,6.},{12.,12.},"I",{8,32},{gus_dir()+"tycho_rave_nside8_Iband.fits",gus_dir()+"tycho_rave_nside32_Iband.fits"}){}
};

// #include "utils.h"
// #include "healpix_base.h"
// #include "healpix_map.h"
// #include "healpix_map_fitsio.h"
// #include "fitshandle.h"

// class RAVE_selection_function{
// private:
//     const int N_I = 70;
//     const int NSIDE = 32;
//     const std::string filename = "RAVE_sf/RAVE_completeness.fits";
//     VecDoub Imag;
//     std::vector<Healpix_Map<double>> HP_Map_grid;
//     int npix;
// public:
//     RAVE_selection_function();
//     int get_number_of_pixels(void){return npix;}
//     double evaluate(double l, double b, double I);
//     double evaluate_at_pixel(double I, int pix, double* l, double* b);
// };


#endif
