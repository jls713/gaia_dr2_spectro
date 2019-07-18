#ifndef IMF_H
#define IMF_H
//=============================================================================
#include <map>
#include <string>
#include <cmath>
#include <memory>
//=============================================================================
#include "params.h"
#include "GSLInterface/GSLInterface.h"
//=============================================================================
/**
 * @brief Initial Mass Function base class
 */
class InitialMassFunction{
protected:
	double MinimumMass,MaximumMass; // minimum and maximum possible stellar masses
	double Norm;
public:
	InitialMassFunction(double MinimumMass,double MaximumMass)
	   :MinimumMass(MinimumMass),MaximumMass(MaximumMass){}
	InitialMassFunction(ModelParameters M);
	/**
	 * @brief initial mass function
	 * @details evaluates initial mass function
	 *
	 * @param m mass in solar masses
	 * @return imf
	 */
	virtual double operator()(double m)=0;
	/**
	 * @brief normalization
	 * @details evaluates \int_{m_min}^{m_max} imf(m) dm
	 * @return normalization
	 */
	double norm(void);
	/**
	 * @brief mean mass
	 * @details evaluates \int_{m_min}^{m_max} m imf(m) dm
	 * @return mean mass
	 */
	double mean_mass(void);
	double sample(void);
};
//=============================================================================

class SalpeterIMF: public InitialMassFunction{
/** Salpeter (1955) IMF **/
public:
	SalpeterIMF(double MinimumMass=0.5, double MaximumMass=100.)
		:InitialMassFunction(MinimumMass,MaximumMass){Norm=1.;Norm=mean_mass();}
	SalpeterIMF(ModelParameters M):InitialMassFunction(M){Norm=1.;Norm=mean_mass();}
	double operator()(double m);
};
class TinsleyIMF: public InitialMassFunction{
/** Tinsley (1980) IMF **/
public:
	TinsleyIMF(double MinimumMass=0.5, double MaximumMass=100.)
		:InitialMassFunction(MinimumMass,MaximumMass){Norm=1.;Norm=mean_mass();}
	TinsleyIMF(ModelParameters M):InitialMassFunction(M){Norm=1.;Norm=mean_mass();}
	double operator()(double m);
};
class ScaloIMF: public InitialMassFunction{
/** Scalo (1998) IMF **/
public:
	ScaloIMF(double MinimumMass=0.5, double MaximumMass=100.)
		:InitialMassFunction(MinimumMass,MaximumMass){Norm=1.;Norm=mean_mass();}
	ScaloIMF(ModelParameters M):InitialMassFunction(M){Norm=1.;Norm=mean_mass();}
	double operator()(double m);
};
class KroupaIMF: public InitialMassFunction{
/** Kroupa (2001) IMF **/
public:
	KroupaIMF(double MinimumMass=0.5, double MaximumMass=100.)
		:InitialMassFunction(MinimumMass,MaximumMass){Norm=1.;Norm=mean_mass();}
	KroupaIMF(ModelParameters M):InitialMassFunction(M){Norm=1.;Norm=mean_mass();}
	double operator()(double m);
};
class KroupaToutGilmoreIMF: public InitialMassFunction{
/** Kroupa, Tout & Gilmore (1993) IMF **/
public:
	KroupaToutGilmoreIMF(double MinimumMass=0.5, double MaximumMass=100.)
		:InitialMassFunction(MinimumMass,MaximumMass){Norm=1.;Norm=mean_mass();}
	KroupaToutGilmoreIMF(ModelParameters M):InitialMassFunction(M){Norm=1.;Norm=mean_mass();}
	double operator()(double m);
};
class ChabrierIMF: public InitialMassFunction{
/** Chabrier (2003) IMF **/
public:
	ChabrierIMF(double MinimumMass=0.5, double MaximumMass=100.)
		:InitialMassFunction(MinimumMass,MaximumMass){Norm=1.;Norm=mean_mass();}
	ChabrierIMF(ModelParameters M):InitialMassFunction(M){Norm=1.;Norm=mean_mass();}
	double operator()(double m);
};
//=============================================================================
// Map for creating new instances of IMF from a string giving the class name
// imf_types[std::string s="<ClassName>"](ModelParameters M)
// produces a shared pointer of Class <ClassName> initialized with the model
// parameters M.
extern shared_map<InitialMassFunction,ModelParameters> imf_types;
//=============================================================================
#endif
//=============================================================================
