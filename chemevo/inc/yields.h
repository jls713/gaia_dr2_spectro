#ifndef YIELDS_H
#define YIELDS_H
//=============================================================================
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <map>
#include <sstream>
#include <memory>
#include <typeinfo>
//=============================================================================
#include "easylogging++.h"
//=============================================================================
#include "utils.h"
#include "params.h"
#include "solar.h"
//=============================================================================
/**
 * @brief remnant mass of initial mass initial_m
 *
 * @param initial_m initial mass
 * @return remnant mass
 */
double remnant_mass(double initial_m);

const double Mch = 1.374999593719; // Chandrasekhar mass
//=============================================================================
/**
 * @brief Base class for storing yields
 */
class Yields{
protected:
	VecDoub yields;		// yields = m_fresh/initial_mass
	VecDoub mass_eject;	// total element ejected mass
	double m_remnant;   // mass of remnant
	double m_eject;		// total mass ejected
public:
	inline double mass(Element E){return mass_eject[element_index[E]];}
	inline double yield(Element E){return yields[element_index[E]];}
	inline double mass_remnant(void){return m_remnant;}
	inline double mass_ejected(void){return m_eject;}
};

//=============================================================================
/**
 * @brief Type Ia Yields
 * @details ejected masses for Type Ia supernovae. No yields given.
 */
class TypeIaYields: public Yields{
private:
	const std::string data_file = "yields/type_Ia/maeda2010.dat";
public:
	TypeIaYields(ModelParameters M);
	inline double yield(Element E) = delete;
};

//=============================================================================
//=============================================================================
/**
 * @brief Yields for a star of mass M and metallicity Z that will form Type II
 */
class TypeIIYields_Single: public Yields{
private:
	double M, Z; // mass and metallicity
public:
	TypeIIYields_Single(VecDoub mass_eject,VecDoub yields,double z, double m,double total_mass_eject,double mass_remnant);
	inline double mass_star(void){return M;}
	inline double metallicity(void){return Z;}
};

//=============================================================================
/**
 * @brief Grid of yields of Type II
 * @details set of TypeIIYields_Single for range of mass and metallicities
 */
class TypeIIYields{
protected:
	VecDoub mass_star, metal;
	std::vector<TypeIIYields_Single> yields;
	AndersSolarAbundances solar;
	double interp_switch(Element E,int index,std::string s);
	/**
	 * @brief general interpolation function
	 * @details interpolates either yields (m_or_y=="yields"), mass for
	 * element E, remnant mass ("mass_remnant") or ejected mass
	 * ("mass_ejected")
	 */
	double interp(Element E, double Mass, double Z,std::string m_or_y="mass");
public:
	TypeIIYields(ModelParameters M):solar(M){}
	/**
	 * @brief ejected mass for element E for star of mass M and metallicity Z
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	inline double mass(Element E,double Mass, double Z){
		return interp(E,Mass,Z,"mass");
	}
	/**
	 * @brief yield for element E for star of mass M and metallicity Z
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	inline double yield(Element E,double Mass, double Z){
		return interp(E,Mass,Z,"yield");
	}
	/**
	 * @brief mass of gas ejected for star of mass M and metallicity Z
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	inline double mass_ejected(double Mass, double Z){
		return interp("",Mass,Z,"mass_ejected");
	}
	/**
	 * @brief mass of remnant of star of mass M and metallicity Z
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	inline double mass_remnant(double Mass, double Z){
		return interp("",Mass,Z,"mass_remnant");
	}
};
/**
 * @brief Kobayashi 2006 yields
 */
class TypeIIKobayashi2006: public TypeIIYields{
private:
	const std::string data_file = "yields/type_II/kobayashi2006.dat";
	const double small_metal = 1e-6; // non-zero metallicity
public:
	TypeIIKobayashi2006(ModelParameters M);
};
/**
 * @brief Chieffi & Limongi 2004 yields
 */
class TypeIIChieffiLimongi2004: public TypeIIYields{
private:
	const std::string data_file = "yields/type_II/chieffi_limongi.dat";
	const double small_metal = 1e-8; // non-zero metallicity
public:
	TypeIIChieffiLimongi2004(ModelParameters M);
};
//=============================================================================
//=============================================================================
/** SUPER AGB YIELDS CURRENTLY NOT IMPLEMENTED. WE ONLY USE TWO SETS OF YIELDS
	AGB AND TYPE II **/
class SuperAGBYields_Single: public Yields{
private:
public:
	SuperAGBYields_Single(std::string yields_string,SolarAbundances &solar);
};

//=============================================================================

class SuperAGBYields{
protected:
	VecDoub mass_star, metal;
	std::vector<SuperAGBYields_Single> yields;
	AndersSolarAbundances solar;
public:
	SuperAGBYields(ModelParameters M):solar(M){}
	virtual double mass(Element E,double Mass, double Z)=0;
	virtual double yield(Element E,double Mass, double Z)=0;
};

//=============================================================================

class SuperAGBYields_None:public SuperAGBYields{
private:
public:
	SuperAGBYields_None(ModelParameters M):SuperAGBYields(M){}
	double mass(Element E,double Mass, double Z){return -1.;}
	double yield(Element E,double Mass, double Z){return -1.;}
};

// class SuperAGBYields_Siess:public SuperAGBYields{
// private:
// 	const std::string data_file = "super_agb/siess2007.dat";
// 	const double small_metal=1e-4;
// public:
// 	SuperAGBYields_Siess(void);
// 	double mass(Element E,double Mass, double Z);
// };

//=============================================================================
//=============================================================================
/**
 * @brief Yields for AGB star of fixed mass and metallicity
 */
class AGBYields_Single: public Yields{
private:
	double M, Z; // mass and metallicity
public:
	AGBYields_Single(void){}
	AGBYields_Single(VecDoub mass_eject,VecDoub yields,double z, double m,double mass_remnant);
	inline double mass_star(void){return M;}
	inline double metallicity(void){return Z;}
};

//=============================================================================
/**
 * @brief Set of AGB yields
 * @details vectors of AGBYields_Single
 *
 */
class AGBYields{
protected:
	std::vector<VecDoub> mass_star;
	VecDoub metal;
	std::vector<std::vector<AGBYields_Single>> yields;
	AndersSolarAbundances solar;
public:
	AGBYields(ModelParameters M):solar(M){}
	virtual inline double mass(Element E,double Mass, double Z){return 0.;}
	/**
	 * @brief yield for element E for star of mass M and metallicity Z
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double yield(Element E,double Mass, double Z){return 0.;}
	/**
	 * @brief mass of gas ejected for star of mass M and metallicity Z
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double mass_ejected(double Mass, double Z){return 0.;}
	/**
	 * @brief mass of remnant of star of mass M and metallicity Z
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double mass_remnant(double Mass, double Z){return 0.;}
};
class AGBYieldsKarakas:public AGBYields{
private:
	const std::string data_file = "yields/agb/karakas2010.dat";

	double interp_switch(Element E, unsigned iz,unsigned im, std::string s);
	/**
	 * @brief general interpolation function
	 * @details interpolates either yields (m_or_y=="yields"), mass for
	 * element E, total ejected mass ("mass_ejected") or remnant mass
	 * ("mass_remnant")
	 */
	double interp(Element E, double Mass, double Z,std::string m_or_y="mass");
public:
	AGBYieldsKarakas(ModelParameters M);
	/**
	 * @brief mass for element E for star of mass M and metallicity Z
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	inline double mass(Element E,double Mass, double Z){
		return interp(E,Mass,Z,"mass");
	}

	/**
	 * @brief yield for element E for star of mass M and metallicity Z
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	inline double yield(Element E,double Mass, double Z){
		return interp(E,Mass,Z,"yield");
	}
	/**
	 * @brief mass of gas ejected for star of mass M and metallicity Z
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	inline double mass_ejected(double Mass, double Z){
		return interp("",Mass,Z,"mass_ejected");
	}
	/**
	 * @brief mass of remnant of star of mass M and metallicity Z
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	inline double mass_remnant(double Mass, double Z){
		return interp("",Mass,Z,"mass_remnant");
	}
};

//=============================================================================
//=============================================================================
/**
 * @brief Full set of yields
 * @details collection of AGB, Super AGB, Type II and Type Ia yields
 * Generic interface class
 * NOTE: SUPER AGB NOT CURRENTLY CONSIDERED
 */
class YieldsSet{
private:
	ModelParameters M;
	std::unique_ptr<AGBYields> agb;
	std::unique_ptr<SuperAGBYields> super_agb;
	std::unique_ptr<TypeIaYields> typeIa;
	std::unique_ptr<TypeIIYields> typeII;
public:
	YieldsSet(ModelParameters M);
	/**
	 * @brief yield for element E for star of mass M and metallicity Z
	 * 		  Includes both Type II and AGB
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	double yield(Element E, double Mass, double Z);
	/**
	 * @brief ejected mass for element E for star of mass M and metallicity Z
	 * 		  Includes both Type II and AGB
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	double mass(Element E, double Mass, double Z);
	/**
	 * @brief remnant mass of star of mass M and metallicity Z
	 * 		  Includes both Type II and AGB
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	double mass_remnant(double Mass, double Z);
	/**
	 * @brief total ejected mass for star of mass M and metallicity Z
	 * 		  Includes both Type II and AGB
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	double mass_ejected(double Mass, double Z);
	/**
	 * @brief Ejected mass for Type Ia SN
	 *
	 * @param E element
	 */
	double typeIa_ejectedmass(Element E);
};

//=============================================================================
// Create pointer to AGB class
template<typename T> AGBYields * createAGBInstance(ModelParameters M) {
	return new T(M);}
// Create pointer to Super AGB class
template<typename T> SuperAGBYields * createSAGBInstance(ModelParameters M) {
	return new T(M);}
// Create pointer to Type II class
template<typename T> TypeIIYields * createTypeIIInstance(ModelParameters M) {
	return new T(M);}
// Create pointer to Type Ia class
template<typename T> TypeIaYields * createTypeIaInstance(ModelParameters M) {
	return new T(M);}
// Map for creating shared pointer instances of AGB from string of class name
extern std::map<std::string, AGBYields*(*)(ModelParameters)> agb_types;
// Map for creating shared pointer instances of SAGB from string of class name
extern std::map<std::string, SuperAGBYields*(*)(ModelParameters)> sagb_types;
// Map for creating shared pointer instances of TypeII from string of class name
extern std::map<std::string, TypeIIYields*(*)(ModelParameters)> typeII_types;
// Map for creating shared pointer instances of TypeIa from string of class name
extern std::map<std::string, TypeIaYields*(*)(ModelParameters)> typeIa_types;
//=============================================================================
#endif
//=============================================================================
