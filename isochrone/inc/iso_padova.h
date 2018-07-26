#ifndef ISOPADOVA_H
#define ISOPADOVA_H
#include "iso_base.h"
#include <map>
//=============================================================================
class isochrone_padova: public isochrone{
private:
    MagMaps mag_maps;
public:
    isochrone_padova(void);
   double get_metallicity(std::vector<std::string> input_iso, std::string dir);
    void fill(std::vector<std::string> input_iso, std::string dir,
              double Age=0., double thin_mag=-1.);
};
//=============================================================================
#endif
//=============================================================================
