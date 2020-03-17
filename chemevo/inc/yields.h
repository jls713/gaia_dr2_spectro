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
#include <regex>
//=============================================================================
#include "easylogging++.h"
//=============================================================================
#include "utils.h"
#include "params.h"
#include "solar.h"
//=============================================================================
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
	void rescale(SolarAbundances *input,
	             std::shared_ptr<SolarAbundances> output){
		for(unsigned i=0;i<mass_eject.size();++i)
			mass_eject[i]+=(output->mass_fraction(i)-input->mass_fraction(i))*mass_ejected();
	}
};

//=============================================================================
/**
 * @brief Type Ia Yields
 * @details ejected masses for Type Ia supernovae. No yields given.
 */
class TypeIaYields: public Yields{
public:
	TypeIaYields(ModelParameters M);
	inline double yield(Element E) = delete;
};
//=============================================================================
/**
 * @brief Maeda 2010 Type Ia Yields
 * @details ejected masses for Type Ia supernovae. No yields given.
 */
class Maeda_TypeIaYields: public TypeIaYields{
private:
	const std::string data_file = "yields/type_Ia/maeda2010.dat";
public:
	Maeda_TypeIaYields(ModelParameters M);
};

//=============================================================================
/**
 * @brief Iwamoto 1999 Type Ia Yields
 * @details ejected masses for Type Ia supernovae. No yields given.
 */
class Iwamoto_TypeIaYields: public TypeIaYields{
private:
	const std::string data_file = "yields/type_Ia/iwamoto_1999.dat";
	const std::string use_model = "W7";
public:
	Iwamoto_TypeIaYields(ModelParameters M);
};

//=============================================================================
/**
 * @brief Seitenzahl 2013 Type Ia Yields
 * @details ejected masses for Type Ia supernovae. No yields given.
 */
class Seitenzahl_TypeIaYields: public TypeIaYields{
private:
	const std::string data_file = "yields/type_Ia/seitenzahl_2013.dat";
public:
	Seitenzahl_TypeIaYields(ModelParameters M);
};
//=============================================================================
/**
 * @brief Thielemann 2003 Type Ia Yields
 * @details ejected masses for Type Ia supernovae. No yields given.
 */
class Thielemann_TypeIaYields: public TypeIaYields{
private:
	const std::string data_file = "yields/type_Ia/thielemann_2003.dat";
public:
	Thielemann_TypeIaYields(ModelParameters M);
};
//=============================================================================
//=============================================================================
/**
 * @brief Yields for a star of mass M and metallicity Z (can use for both Type II and AGB)
 */
class Yields_Single: public Yields{
private:
	double M, Z; // mass and metallicity
public:
	Yields_Single(VecDoub mass_eject,VecDoub yields,double z, double m,
	              double total_mass_eject,double mass_remnant);
	inline double mass_star(void){return M;}
	inline double metallicity(void){return Z;}
};
template<class T>
bool comp_function_pair( const std::pair<unsigned, std::pair<T,T>>& i, const std::pair<unsigned,std::pair<T,T>>& j ) {
    if( j.second.first < i.second.first ) return false;
    if( i.second.first < j.second.first ) return true;
    return i.second.second < j.second.second;
}
template<class T>
bool comp_function(const std::pair<unsigned, T>& i,
                   const std::pair<unsigned, T>& j ) {
    return i.second < j.second;
}
template< class T >
void reorder(std::vector<T> &v, std::vector<size_t> const &order )  {
    for ( int s = 1, d; s < order.size(); ++ s ) {
        for ( d = order[s]; d < s; d = order[d] ) ;
        if ( d == s ) while ( d = order[d], d != s ) std::swap( v[s], v[d] );
    }
}
//=============================================================================
/**
 * @brief Regular Grid of yields
 * @details set of Yields_Single for range of mass and metallicities
 */
class Yields_Grid{
// protected:
// 	VecDoub mass_star, metal;
// 	std::vector<Yields_Single> yields;

protected:
	VecDoub mass_star;
	VecDoub metal;
	std::vector<Yields_Single> yields;
	virtual double interp_switch(Element E,int index,std::string s);
	/**
	 * @brief general interpolation function
	 * @details interpolates either yields (m_or_y=="yields"), mass for
	 * element E, remnant mass ("mass_remnant") or ejected mass
	 * ("mass_ejected")
	 */
	virtual double interp(Element E, double Mass, double Z,std::string m_or_y="mass");
    virtual void sort_arrays(void);
	void rescale(SolarAbundances *input,
	             std::shared_ptr<SolarAbundances> output){
		for(auto i=0;i<yields.size();++i)
			yields[i].rescale(input, output);
	}
public:
	Yields_Grid(ModelParameters M){}
	/**
	 * @brief ejected mass for element E for star of mass M and metallicity Z
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double mass(Element E,double Mass, double Z){
		return interp(E,Mass,Z,"mass");
	}
	/**
	 * @brief yield for element E for star of mass M and metallicity Z
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double yield(Element E,double Mass, double Z){
		return interp(E,Mass,Z,"yield");
	}
	/**
	 * @brief mass of gas ejected for star of mass M and metallicity Z
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double mass_ejected(double Mass, double Z){
		return interp("",Mass,Z,"mass_ejected");
	}
	/**
	 * @brief mass of remnant of star of mass M and metallicity Z
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double mass_remnant(double Mass, double Z){
		return interp("",Mass,Z,"mass_remnant");
	}
};
class Yields_Grid_None:public Yields_Grid{
public:
	Yields_Grid_None(ModelParameters M): Yields_Grid(M){}
	inline double mass(Element E,double Mass, double Z){ return 0.;}
	inline double yield(Element E,double Mass, double Z){ return 0.;}
	inline double mass_ejected(double Mass, double Z){ return 0.;}
	inline double mass_remnant(double Mass, double Z){ return 0.;}
};
class TypeIIYields: public Yields_Grid{
protected:
	TypeIIYields(ModelParameters M): Yields_Grid(M){}
};
/**
 * @brief Kobayashi 2006 yields
 */
class TypeIIKobayashi2006: public TypeIIYields{
private:
	const std::string data_file = "yields/type_II/kobayashi2006.dat";
	const double small_metal = 1e-6; // non-zero metallicity
	AndersSolarAbundances solar;
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
	AndersSolarAbundances solar;
public:
	TypeIIChieffiLimongi2004(ModelParameters M);
};
/**
 * @brief Chieffi & Limongi 2018 yields
 */
class TypeIIChieffiLimongi2018: public TypeIIYields{
private:
	const std::string data_file = "yields/type_II/chieffi_limongi2018.dat";
	const int vel_target = 0;
	std::map<double, double> metal_list = {{0,1.345e-2},{-1,3.236e-3},
												{-2,3.236e-4},{-3,3.236e-5}};
	VecDoub mass_list = {13.,15.,20.,25.,30.,40.,60.,80.,120.};
	AsplundSolarAbundances solar;
public:
	TypeIIChieffiLimongi2018(ModelParameters M);
};
/**
 * @brief Nugrid yields from Pignatari 2016 and Ritter 2018
 */
class NugridYields: public Yields_Grid{
private:
	const std::string data_file = "yields/nugrid_set1ext.dat";
	GrevesseNoelSolarAbundances solar;
public:
	NugridYields(ModelParameters M);
};
//=============================================================================
/**
 * @brief Set of yields but for an irregular gridding
 * @details Set of metallicity-dependent yields but each metallicity has a
 * 			different range of masses. vectors of Yields_Single
 *
 */
class Yields_Grid_Irregular: public Yields_Grid{
protected:
    std::vector<VecDoub> mass_star_vec;
    std::vector<std::vector<Yields_Single>> yields_vec;
	void sort_arrays(void);

	double interp_switch(Element E, unsigned iz,unsigned im, std::string s);
	/**
	 * @brief general interpolation function
	 * @details interpolates either yields (m_or_y=="yields"), mass for
	 * element E, total ejected mass ("mass_ejected") or remnant mass
	 * ("mass_remnant")
	 */
	double interp(Element E, double Mass, double Z,std::string m_or_y="mass");
	void rescale(SolarAbundances *input,
	             std::shared_ptr<SolarAbundances> output){
		for(auto i=0;i<mass_star_vec.size();++i)
			for(auto j=0;j<mass_star_vec[i].size();++j)
				yields_vec[i][j].rescale(input, output);
	}
public:
	Yields_Grid_Irregular(ModelParameters M):Yields_Grid(M){}
	/**
	 * @brief mass for element E for star of mass M and metallicity Z
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double mass(Element E,double Mass, double Z){
		return interp(E,Mass,Z,"mass");
	}

	/**
	 * @brief yield for element E for star of mass M and metallicity Z
	 *
	 * @param E element
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double yield(Element E,double Mass, double Z){
		return interp(E,Mass,Z,"yield");
	}
	/**
	 * @brief mass of gas ejected for star of mass M and metallicity Z
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double mass_ejected(double Mass, double Z){
		return interp("",Mass,Z,"mass_ejected");
	}
	/**
	 * @brief mass of remnant of star of mass M and metallicity Z
	 *
	 * @param Mass mass
	 * @param Z metallicity
	 */
	virtual inline double mass_remnant(double Mass, double Z){
		return interp("",Mass,Z,"mass_remnant");
	}
};
class AGBYieldsKarakas:public Yields_Grid_Irregular{
public:
	AndersSolarAbundances solar;
private:
	const std::string data_file = "yields/agb/karakas2010.dat";
public:
	AGBYieldsKarakas(ModelParameters M);
};
class AGBYieldsVentura:public Yields_Grid_Irregular{
public:
	GrevesseSauvalSolarAbundances solar;
private:
	const std::string data_file = "yields/agb/Ventura_2013.dat";
public:
	AGBYieldsVentura(ModelParameters M);
};
//=============================================================================
//=============================================================================
/**
 * @brief Full set of yields
 * @details collection of AGB, Type II and Type Ia yields
 * Generic interface class
 */
class YieldsSet{
private:
	ModelParameters M;
	std::unique_ptr<Yields_Grid> agb;
	std::unique_ptr<TypeIaYields> typeIa;
	std::unique_ptr<Yields_Grid> typeII;
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
template<typename T> Yields_Grid * createAGBInstance(ModelParameters M) {
	return new T(M);}
// Create pointer to Type II class
template<typename T> Yields_Grid * createTypeIIInstance(ModelParameters M) {
	return new T(M);}
// Create pointer to Type Ia class
template<typename T> TypeIaYields * createTypeIaInstance(ModelParameters M) {
	return new T(M);}

// Map for creating shared pointer instances of AGB from string of class name
extern std::map<std::string, Yields_Grid*(*)(ModelParameters)> agb_types;
// Map for creating shared pointer instances of TypeII from string of class name
extern std::map<std::string, Yields_Grid*(*)(ModelParameters)> typeII_types;
// Map for creating shared pointer instances of TypeIa from string of class name
extern std::map<std::string, TypeIaYields*(*)(ModelParameters)> typeIa_types;
//=============================================================================
#endif
//=============================================================================
