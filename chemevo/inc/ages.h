#ifndef AGES_H
#define AGES_H
//=============================================================================
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <memory>
//=============================================================================
#include "utils.h"
#include "GSLInterface/GSLInterface.h"
#include "params.h"
#include "easylogging++.h"
//=============================================================================
/**
 * @brief Base class for storing stellar lifetime functions
 */
class StellarLifetime{
public:
	/**
	 * @brief stellar lifetime as a function of mass m and metallicity Z
	 *
	 * @param m mass
	 * @param Z metallicity
	 *
	 * @return lifetime
	 */
	virtual double operator()(double m, double Z)=0;
	/**
	 * @brief mass of a star that has age t and metallicity Z
	 *
	 * @param t age
	 * @param Z metallicity
	 *
	 * @return mass of star
	 */
	virtual double mass_star_dying_now(double t, double Z)=0;
	/**
	 * @brief gradient of log M wrt log t for stars of age t and metallicity Z
	 *
	 * @param t age
	 * @param Z metallicity
	 *
	 * @return gradient of log M wrt log t
	 */
	virtual double dlogMdlogt(double t, double Z)=0;
};
//=============================================================================
class MaederMeynet1989:public StellarLifetime{
public:
	MaederMeynet1989(ModelParameters M){}
	double operator()(double m, double Z);
	double mass_star_dying_now(double t, double Z);
	double dlogMdlogt(double t, double Z);
};

// class Tinsley1980{
// public:
// 	Tinsley1980(void){}
// 	double operator()(double m, double Z);
// };

class Kodama1997:public StellarLifetime{
public:
	Kodama1997(ModelParameters M){}
	double operator()(double m, double Z);
	double mass_star_dying_now(double t, double Z);
	double dlogMdlogt(double t, double Z);
};

class PadovaniMatteucci1993:public StellarLifetime{
public:
	PadovaniMatteucci1993(ModelParameters M){}
	double operator()(double m, double Z);
	double mass_star_dying_now(double t, double Z);
	double dlogMdlogt(double t, double Z);
};

class InterpolateAgeGrid: public StellarLifetime{
protected:
	const VecDoub metal_grid;// log(Z)
	VecDoub mass_grid, age_grid;
	std::vector<VecDoub> lifetimes, masses;
	std::shared_ptr<interpolator2D> int2D, int2D_back;
public:
	InterpolateAgeGrid(VecDoub mg): metal_grid(mg){}
	double operator()(double m, double Z);
	double mass_star_dying_now(double t, double Z);
	double mass_star_dying_now2(double t, double Z);
	double dlogMdlogt(double t, double Z);
};

class Portinari1998:public InterpolateAgeGrid{
protected:
	const std::string data_file = "portinari_ages.dat";
	// const VecDoub metal_grid
	// 		= {-3.398,-2.398,-2.0969,-1.69897,-1.3010};// log(Z)
public:
	Portinari1998(ModelParameters M); // Load grid
};

class PadovaIsochrones:public InterpolateAgeGrid{
protected:
	const std::string data_file = "padova_ages.dat";
	// const VecDoub metal_grid
	// 		= create_log_range(-3., 0.7, 20);// log(Z)
public:
	PadovaIsochrones(ModelParameters M); // Load grid
};
//=============================================================================
// Map for creating new instances of stellar lifetimes from a string giving
// the class name life_types[std::string s="<ClassName>"](ModelParameters M)
// produces a shared pointer of Class <ClassName> initialized with the model
// parameters M.
extern shared_map<StellarLifetime,ModelParameters> life_types;
//=============================================================================
#endif
//=============================================================================
