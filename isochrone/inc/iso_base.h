#ifndef ISOBASE_H
#define ISOBASE_H
//=============================================================================
#include "utils_iso.h"
#include "extinction.h"
#include <map>
#include <utility>
#include <vector>

//=============================================================================
namespace SUNDATA{
    // Data on the Sun, metallicity Z, Helium abundance Y,
    // effective temperature and log g
    const double SUNZ = .0198;
    const double SUNY= .273;
    const double SUNTEFF = 3.75935;
    const double SUNLOGG = 4.44;
}
//=============================================================================
typedef std::map<std::string, VecDoub> MagData;
typedef std::map<std::string,
                 std::vector<std::pair<std::string, unsigned> > > MagMaps;
//=============================================================================
/**
 * @brief Base class for storing a single isochrone i.e. a track of initial
 * mass against various stellar properties for fixed metallicity and age
 */
class isochrone{
    protected:
        int N_length; 				// number of isochrone points
        double FeH, age;			// metallicity and age of isochrone
        VecDoub InitialMass;		// Grid of the initial mass
        VecDoub Mass, Teff, L, Logg;// Grids of properties
                                    // (mass, effective temperature,
                                    //  luminosity and surface gravity)
        VecDoub es;                 // = 1-age/maxage
        MagData mags;               // stores magnitude data
        VecDoub mf;                 // Mass/MaxMass
    public:
		//=====================================================================
		// Simple print
        inline void print(){
        	std::cout<<"FeH = "<<FeH<<", age = "<<age<<std::endl;}
		//=====================================================================
		// Getters
        inline int N(){
        	return N_length;
        }
        inline double feh(){
        	return FeH;
        }
        inline double tau(){
        	return age;
        }
        inline double initial_mass(int j){
        	return InitialMass[j];
        }
        VecDoub const &im() const {
        	return InitialMass;
        }
        inline double mass(int j){
        	return Mass[j];
        }
        inline double logTeff(int j){
        	return Teff[j];
        }
        inline double logg(int j){
        	return Logg[j];
        }
        inline double logg_calc(double M, double LL, double TT){
			return log10(M)-LL+4.*(TT-  SUNDATA::SUNTEFF)+SUNDATA::SUNLOGG;
        }
        inline double maxmass(){
        	return InitialMass.back();
        }
        inline double minmass(){
        	return InitialMass.front();
        }
        VecDoub const &afraction() const {
            return es;
        }
        inline double age_fraction(int j){
            return es[j];
        }
        inline void set_maxage(unsigned m,double a){
            es[m]=age/a;
        }
        /**
         * @brief initial mass separation between two adjacent isochrone points
         * @details distance between isochrone point index_m and its two
         * nearest neighbours
         *
         * @param index_m isochrone point index
         * @return \Delta M
         */
        inline double delta_mass(int index_m){
		    // find the volume occupied by an initial mass point
		    if(index_m==0)return fabs(InitialMass[1]-InitialMass[0]);
		    if(index_m==N_length-1)return fabs(InitialMass[N_length-1]-InitialMass[N_length-2]);
		    else return fabs((InitialMass[index_m+1]-InitialMass[index_m-1])/2.);
		}
        inline double mag(int j, std::string mag){
            return mags[mag][j];
        }
        inline double colour(int j, std::vector<std::string> colour){
            return mags[colour[0]][j]-mags[colour[1]][j];
        }
        inline double distance(int j, double M, std::string mag){
            return pow(10.,0.2*(mags[mag][j]-M)-2.);
        }
};

//=============================================================================
#endif
//=============================================================================
