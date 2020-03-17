#ifndef SOLAR_H
#define SOLAR_H
//=============================================================================
#include "utils.h"
#include "params.h"
//=============================================================================
typedef std::string Element;
extern std::map<Element,int> element_index; // mapping element to atomic number
extern std::map<int,Element> element_index_reverse; // reverse of above
extern std::map<Element,bool> is_alpha_element;// returns true if alpha element
//=============================================================================
/**
 * @brief Solar abundance data
 * @details storage class for solar abundances
 *
 */
class SolarAbundances{
protected:
	const double X_p=0.7547; // Primordial (cosmological) Hydrogen fraction
	const double Y_p=0.2453; // Primordial (cosmological) Helium fraction
	VecDoub masses;			// Element masses in sun
	VecDoub abundances;		// Abundances in sun
	VecDoub mfracs;			// Mass fractions in sun
	double X_, Y_, Z_;		// Hydrogen, Helium and Metals mass fraction
public:

	// SolarAbundances(ModelParameters M){}

	double mass_fraction(unsigned e){
		return mfracs[e];
	}
	double mass_fraction(Element E){
		return mass_fraction(element_index[E]);
	}
	double abundance(unsigned e){
		return abundances[e];
	}
	double abundance(Element E){
		return abundance(element_index[E]);
	}
	/**
	 * @brief Scaled solar mass fraction
	 * @details scale mass fractions relative to solar using (Z/Z_sun). Note that for X and Y the zero-points are X_p and Y_p (the primordial values).
	 *
	 * @param E Element
	 * @param z metallicity
	 *
	 * @return scale solar abundance
	 */
    double scaled_solar_mass_frac(Element E, double z, double alpha_enhance=0.){
    	double out=0.;
		if(E=="H") out= X_p-(X_p-X())/Z()*z;
		else if(E=="He") out= Y_p+(Y()-Y_p)/Z()*z;
		else out = mass_fraction(E)*z/Z();
		if(is_alpha_element[E])
			out*=pow(10., alpha_enhance);
		return out;
	}
	/**
	 * @brief Scaled solar mass fraction
	 * @details scale mass fractions relative to solar using (Z/Z_sun). Note that for X and Y the zero-points are X_p and Y_p (the primordial values).
	 *
	 * @param e atomic number
	 * @param z metallicity
	 *
	 * @return scale solar abundance
	 */
    double scaled_solar_mass_frac(unsigned e, double z, double alpha_enhance=0.){
    	double out=0.;
		if(e==0) out = X_p-(X_p-X())/Z()*z;
		else if(e==1) out = Y_p+(Y()-Y_p)/Z()*z;
		else out = mass_fraction(e)*z/Z();
		if(is_alpha_element[element_index_reverse[e]])
			out*=pow(10., alpha_enhance);
		return out;
	}
    inline double X(void){return X_;}
    inline double Y(void){return Y_;}
    inline double Z(void){return Z_;}
};
//=============================================================================
// Abundances from Asplund (2008)
class AsplundSolarAbundances: public SolarAbundances{
private:
	const std::string data_file = "asplund_abundances.dat";
public:
	AsplundSolarAbundances(ModelParameters M);
};
// Generic loader for abundances for Grevesse & Noel (1993) and Anders & Grevesse 1989
class AndersTypeSolarAbundances: public SolarAbundances{
public:
	AndersTypeSolarAbundances(ModelParameters M, std::string data_folder);
};
// Abundances from Anders & Grevesse 1989
class AndersSolarAbundances: public AndersTypeSolarAbundances{
public:
	AndersSolarAbundances(ModelParameters M)
		: AndersTypeSolarAbundances(M, "anders1989_abundances.dat"){}
};
// Abundances from Grevesse & Noel (1993) and Anders & Grevesse 1989
class GrevesseNoelSolarAbundances: public AndersTypeSolarAbundances{
public:
	GrevesseNoelSolarAbundances(ModelParameters M)
		: AndersTypeSolarAbundances(M, "grevesse_noel1993_abundances.dat"){}
};
// Abundances from Grevesse & Sauval (1998) and Anders & Grevesse 1989
class GrevesseSauvalSolarAbundances: public AndersTypeSolarAbundances{
public:
	GrevesseSauvalSolarAbundances(ModelParameters M)
		: AndersTypeSolarAbundances(M, "grevesse_sauval1998abundances.dat"){}
};
//=============================================================================
extern shared_map<SolarAbundances,ModelParameters> solar_types;
//=============================================================================
#endif
//=============================================================================
