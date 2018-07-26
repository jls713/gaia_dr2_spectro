#ifndef ISOGRID_H
#define ISOGRID_H
//=============================================================================
#include "json.hpp"
#include "iso_base.h"
#include "iso_basti.h"
#include "iso_dartmouth.h"
#include "iso_padova.h"
//=============================================================================
/**
 * @brief Class for storing a grid of isochrones
 *
 * @tparam isochrone_g isochrone type for grid
 */
template <class isochrone_g>
class isochrone_grid{
    private:
    public:
        int NF, NA;    // Number of metallicities and ages
        int Np;        // Number of points along isochrone
        unsigned int iso_grid_size;       // total number of isochrones
        std::vector<isochrone_g> iso_grid;// grid for storing isochrones
        VecDoub agegrid, fehgrid;         // ages and metallicities of isos
        std::vector<VecDoub> massmax; // maximum mass at age and metallicity
        //=====================================================================
        /**
         * @brief load in a grid of isochrones
         *
         * @param type -- type of isochrone (Sloan, Johnson, Johnson_More, Padova_dense, Dartmouth)
         * @param thin -- store every thin entry in isochrone
         * @param feh_err -- down sample metallicities until spacing=feh_err/2
         * @param thin_mag -- down sample magnitude spacing to less than thin_mag
         */
        isochrone_grid(std::string type = "Padova",int thin=1, double feh_err=1e-4, double thin_mag=-1.);
        /**
         * @brief separation in age between neighbouring isochrones
         *
         * @param index_a location of isochrone in age
         */
        inline double delta_age(int index_a){
            // find the volume occupied by an age point
            if(index_a==0)return fabs(agegrid[1]-agegrid[0]);
            if(index_a==NA-1)return fabs(agegrid.back()-agegrid[NA-2]);
            else return fabs((agegrid[index_a+1]-agegrid[index_a-1])/2.);
        }
        /**
         * @brief separation in metallicity between neighbouring isochrones
         *
         * @param index_f location in metallicity
         */
        inline double delta_Z(int index_f){
            // find the volume occupied by an Z point
            if(index_f==0)return fabs(fehgrid[1]-fehgrid[0]);
            if(index_f==NF-1)return fabs(fehgrid.back()-fehgrid[NF-2]);
            else return fabs((fehgrid[index_f+1]-fehgrid[index_f-1])/2.);
        }
        /**
         * @brief isochrone at location (i,j) in grid
         *
         * @param i age location
         * @param j metallicity location
         *
         * @return isochrone at (i,j)
         */
        isochrone_g *iso(int i,int j);
        /**
         * @brief indices of isochrone nearest to metallicity, age and mass
         *
         * @param Z metallicity
         * @param age age
         * @param M mass
         * @return indices of nearest isochrone
         */
        std::vector<int> find_nearest(double Z, double age, double M);
        /**
         * @brief properties of isochrone interpolated on grid
         * @detail returns (Teff,logg,mag bands vec) for interpolated isochrone
         *
         * @param Z metallicity
         * @param age age
         * @param M initial mass
         * @param band vector of magnitude bands requried
         */
        VecDoub interp(double Z, double age, double M, std::vector<std::string> band);
        VecDoub interp_es(double Z, double age, double M, std::vector<std::string> band);
        /**
         * @brief maximum mass at given age and metallicity
         *
         * @param Z metallicity
         * @param age
         *
         * @return max mass
         */
        double max_mass(double Z, double age);
        /**
         * @brief maximum mass at given mass and metallicity
         *
         * @param Z metallicity
         * @param M initial mass
         *
         * @return max age
         */
        double max_age(double Z, double M);
        /**
         * @brief finds the points that bracket (Z,age,es)
         *
         * @param Z metallicity
         * @param age age
         * @param es evolutionary stage = tau/tau_max
         * @param bot_Z lower metallicity grid point
         * @param top_Z upper metallicity
         * @param bot_a lower age grid point
         * @param top_a upper age grid point
         * @param bot_m00 lower mass grid point at lower Z, lower age
         * @param top_m00 etc.
         * @param bot_m01 etc.
         * @param top_m01 etc.
         * @param bot_m10 etc.
         * @param top_m10 etc.
         * @param bot_m11 etc.
         * @param top_m11 etc.
         * @param dZ (Z-Z[lower])/(Z[upper]-Z[lower])
         * @param mdZ 1-dZ
         * @param da etc.
         * @param mda etc.
         * @param dm00 etc.
         * @param mdm00 etc.
         * @param dm10 etc.
         * @param mdm10 etc.
         * @param dm01 etc.
         * @param mdm01 etc.
         * @param dm11 etc.
         * @param mdm11 etc.
         */
        void find_es_intercepts(double Z, double age, double es, int *bot_Z, int *top_Z, int *bot_a, int *top_a, int *bot_m00,int *top_m00,int *bot_m01,int *top_m01,int *bot_m10,int *top_m10,int *bot_m11,int *top_m11,double *dZ,double *mdZ,double *da,double *mda,double *dm00,double *mdm00,double *dm10,double *mdm10,double *dm01,double *mdm01,double *dm11,double *mdm11);
        /**
         * @brief finds initial mass from evolutionary stage = tau/tau_max
         *
         * @param Z metallicity
         * @param age age
         * @param es evo stage
         * @return initial mass
         */
        double mass_from_es(double Z, double age, double es);
        /**
         * @brief finds gradient of evo stage with initial mass.
         *
         * @param Z metallicity
         * @param age age
         * @param es evo stage
         * @return d(es)/dM
         */
        double d_es_d_M(double Z, double age, double es);
        /**
         * @brief thin feh grid until mean spacing < feh_err/2
         */
        void thin_feh_grid(double feh_err);
    };


//=============================================================================
#endif
//=============================================================================
