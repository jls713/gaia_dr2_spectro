#ifndef GALAXY_H
#define GALAXY_H
//=============================================================================
#include "utils.h"
#include "utils_ch.h"
#include "params.h"
//=============================================================================
class Grid{
private:
	VecDoub radial_grid; // radial grid points
	VecDoub time_grid;   // time grid points
	VecDoub age_grid;    // age grid points = MaxAge-time_grid
	std::vector<VecDoub> grid; /** grid in radius and time **/
	double MaxAge;
public:
	//=========================================================================
	// Constructors
	Grid(unsigned NR, unsigned NT, double Rmin, double Rmax, double MaxAge);
	Grid(ModelParameters params);
	//=========================================================================
	// Getters
	/**
	 * @brief evaluates grid at radius R and time t
	 * @details linearly interpolates on grid to find properties at
	 * radius R and time t
	 *
	 * @param R radius
	 * @param t time
	 * @param loginterp , if true interpolate using logarithm
	 * @param extrapolate , if true log extrapolate off ends of R grid
	 *
	 * @return grid properties
	 */
	double operator()(double R, double t, bool loginterp=false, bool extrapolate=false);
	/**
	 * @brief evaluates time gradient at radius R and time t
	 * @details linearly interpolates on grid to find time gradient at
	 * radius R and time t
	 *
	 * @param R radius
	 * @param t time
	 *
	 * @return time gradient
	 */
	double t_gradient(double R, double t);
	/**
	 * @brief evaluates radial gradient at radius R and time t
	 * @details linearly interpolates on grid to find radial gradient at
	 * radius R and time t
	 *
	 * @param R radius
	 * @param t time
	 *
	 * @return radial gradient
	 */
	double R_gradient(double R, double t);
	/**
	 * @brief evaluates grid properties at radius R and ***age*** age
	 * @details differs from above as evaluates at age not time
	 *
	 * @param R radius
	 * @param age
	 *
	 * @return grid properties
	 */
	inline double properties(double R, double age)
		{return (*this)(R,MaxAge-age);}
	inline VecDoub grid_radial(void){ return radial_grid;}
	inline VecDoub grid_time(void){ return time_grid;}
	inline VecDoub grid_age(void){ return age_grid;}
	inline double operator()(unsigned r, unsigned t){return grid[r][t];}
	inline VecDoub grid_fixed_t(unsigned t){ return transpose(grid)[t];}
	inline VecDoub grid_fixed_r(unsigned r){ return grid[r];}
	inline double mid_last_grid_r(){
		unsigned nr=radial_grid.size();
		return .5*(radial_grid[nr-2]+radial_grid[nr-1]);
	}
	double log_extrapolate_low(double R, double t);
	double log_extrapolate_low(double R, unsigned t);
	double log_extrapolate_high(double R, double t);
	double log_extrapolate_high(double R, unsigned t);
	double Rdown(unsigned nR);
	double Rup(unsigned nR);
	double annulus_area(int nR);
	double annulus_width(int nR);
	// insert extra time t column before entry nt
	void add_time(unsigned nt, double t);
	//=========================================================================
	// Scaling all grid entries
	inline void scale(double s, unsigned nR, unsigned nt){grid[nR][nt]*=s;}
	inline void log10scale(double solar){
		for(unsigned nR=0;nR<radial_grid.size();++nR)
			for(unsigned nt=0;nt<time_grid.size();++nt)
				grid[nR][nt]=log10(grid[nR][nt]/solar);
	}
	//=========================================================================
	// Setters
	void set_fixed_r_t_const(double c);
	void set_fixed_t_const(double c,unsigned t);
	void set_fixed_t(VecDoub rad,unsigned t);
	void set_fixed_r(VecDoub tim,unsigned r);
	void set(double x, unsigned r, unsigned t);
	//=========================================================================
	// Write HDF5
	void write_hdf5(H5::H5File &fout,std::string name);
	void write_1D_vector(H5::H5File &fout, VecDoub &v, std::string name);
	void write_radial_grid_hdf5(H5::H5File &fout);
	void write_time_grid_hdf5(H5::H5File &fout);
	//=========================================================================
	// Printing
	void print_grid(std::ostream &out);
	void print_grid_fixed_radius(std::ostream &out,double R);
	void print_grid_fixed_time(std::ostream &out,double t);
	void print_grid_now(std::ostream &out);
};
//=============================================================================
#endif
//=============================================================================
