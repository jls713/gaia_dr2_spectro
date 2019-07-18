#ifndef ONSKY_SF_H
#define ONSKY_SF_H

#include <memory>
#include <map>
#include "utils.h"
#include "healpix_base.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "fitshandle.h"

class Selection_function{
public:
    virtual std::vector<std::string> mag_band(void)=0;
    virtual double evaluate(double l, double b,
                            std::vector<double> mag, bool interp=false)=0;
    virtual double min_mag(void)=0;
    virtual double max_mag(void)=0;
};
class NullSelection_Function:public Selection_function{
public:
    std::vector<std::string> mag_band(void){return {};}
    double evaluate(double l, double b, std::vector<double> mag, bool interp=false){
        return 1.;
    }
    double min_mag(void){return 0.;}
    double max_mag(void){return 0.;}
};
class Binned_selection_function:public Selection_function{
    // returns zero outside mag range supplied
private:
    unsigned N_mag; // Number of magnitude bins
    int NSIDE; // Healpix NSIDE
    std::string band; // magnitude band
    std::string filename; // fits file with selection function
    VecDoub mag_edges, mag;    // bin edges and centres in magnitude
    std::vector<Healpix_Map<double>> HP_Map_grid;
    int npix;
public:
    Binned_selection_function(unsigned N_mag, double magmin, double magmax, std::string band, int NSIDE, std::string filename);
    std::vector<std::string> mag_band(void){return {band};}
    int get_number_of_pixels(void){return npix;}
    double evaluate(double l, double b, std::vector<double> I,
                    bool interp=false);
    double evaluate_at_pixel(double I, int pix, double* l, double* b);
    double pixel_to_lb(int pix, double *l, double *b);
    inline double min_mag(void){ return mag_edges.front();}
    inline double max_mag(void){ return mag_edges.back();}
};
// Two selection function maps for low and high magnitude ranges
class Hybrid_Binned_selection_function: public Selection_function{
private:
    double mag_split;
    std::string band;
    std::unique_ptr<Binned_selection_function> lowMap, highMap;
public:
    Hybrid_Binned_selection_function(double mag_split,VecInt N_mag, VecDoub magmin, VecDoub magmax, std::string band, VecInt NSIDE, std::vector<std::string> filenames):
	mag_split(mag_split),band(band),
	lowMap(new Binned_selection_function(N_mag[0],magmin[0],magmax[0],band,NSIDE[0],filenames[0])),
	highMap(new Binned_selection_function(N_mag[1],magmin[1],magmax[1],band,NSIDE[1],filenames[1])){}
    std::vector<std::string> mag_band(void){return {band};}
	double evaluate(double l, double b, std::vector<double> mag, bool interp=false){
		if(mag[0]<mag_split) return lowMap->evaluate(l,b,mag,interp);
		else return highMap->evaluate(l,b,mag,interp);
	}
    int get_number_of_pixels(void){return highMap->get_number_of_pixels();}
    inline double min_mag(void){ return lowMap->min_mag();}
    inline double max_mag(void){ return highMap->max_mag();}
};

class Stacked_Binned_selection_function: public Selection_function{
private:
    std::vector<Selection_function*> sf;
public:
    Stacked_Binned_selection_function(std::vector<Selection_function*> sf):sf(std::move(sf)){}
    // Return Magnitude band for first selection function
    std::vector<std::string> mag_band(void){
        std::vector<std::string> s;
        for(auto SF: sf)
            s.push_back(SF->mag_band()[0]);
        return s;
    }
    double evaluate(double l, double b, std::vector<double> mag, bool interp=false){
        double e = 1.;int n=0;
        for(auto SF: sf){
            e*=SF->evaluate(l,b,{mag[n]},interp);++n;
        }
        return e;
    }
    inline double min_mag(void){ return sf[0]->min_mag();}
    inline double max_mag(void){ return sf[0]->max_mag();}
};

#endif
