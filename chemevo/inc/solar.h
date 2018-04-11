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
    double scaled_solar_mass_frac(Element E, double z){
		if(E=="H") return X_p-(X_p-X())/Z()*z;
		else if(E=="He") return Y_p+(Y()-Y_p)/Z()*z;
		else return mass_fraction(E)*z/Z();
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
    double scaled_solar_mass_frac(unsigned e, double z){
		if(e==0) return X_p-(X_p-X())/Z()*z;
		else if(e==0) return Y_p+(Y()-Y_p)/Z()*z;
		else return mass_fraction(e)*z/Z();
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
// Abundances from Anders & Grevesse 1989
class AndersSolarAbundances: public SolarAbundances{
private:
	const std::string data_file = "anders1989_abundances.dat";
public:
	AndersSolarAbundances(ModelParameters M);
};
//=============================================================================
#endif
//=============================================================================
