#ifndef UTILS_ISO_H
#define UTILS_ISO_H
//=============================================================================
// A few utility functions for isochrone interface
//=============================================================================
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <regex>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fcntl.h>
#include "json.hpp"
using json = nlohmann::json;
#include "utils.h"
//=============================================================================
/**
 * @brief make unique pointer
 *
 * @param args constructor arguments
 * @tparam T class type to make
 * @tparam typename ...Args class types of constructor arguments
 * @return unique pointer
 */
template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}
//=============================================================================
/**
 * @brief get files in a director
 * @details get all files in a directory and return as vector of strings
 *
 * @param out vector of strings of filenames in directory
 * @param directory string with filename
 */
void GetFilesInDirectory(std::vector<std::string> &out, const std::string &directory);
/**
 * @brief Simple Kroupa initial mass function dN/dM
 *
 * @param M mass in solar masses
 * @return dN/dM
 */
double KroupaIMF_default(double M);
/**
 * @brief Approximate I band from 2MASS J & K bands
 *
 * @param J 2MASS band J
 * @param K 2MASS band Ks
 *
 * @return approximate I band
 */
double I2MASS(double J, double K);
/**
 * @brief J, K & H VISTA bands
 *
 * @param J 2MASS band J
 * @param K 2MASS band Ks
 * @param H 2MASS band H
 *
 * @return Vista JKH
 */
VecDoub JKH_vista_from_2MASS(double J, double K, double H);
/**
 * @brief J, K & H 2MASS bands
 *
 * @param J Johnson band J
 * @param K Johnson band Ks
 * @param H Johnson band H
 *
 * @return 2MASS JKH
 */
VecDoub JKH_2MASS(double J, double K, double H);
/**
 * @brief magnitude of vector
 *
 * @param x input vector
 * @return |x|
 */
double magnitude(VecDoub x);
/**
 * @brief V' magnitude -- approximately straightens selection function in V
 * against (B-V)
 *
 * @param B Johnson B
 * @param V Johnson V
 *
 * @return V+0.57(B-V)
 */
double VpMag(double B, double V);
/**
 * @brief cross product
 *
 * @param a input vector 1
 * @param b input vector 2
 *
 * @return cross product a x b
 */
double cross_product2D(const VecDoub &a, const VecDoub &b);

json parameters(void);
/**
 * @brief split string
 *
 * @param a input string
 * @param b string to split on
 *
 * @return vector of substrings
 */
std::vector<std::string> Split(const std::string& subject, const std::string& split_string);

// Custom comparator for Padova isochrone files
// Sorts alphabetical and then numerical
// e.g. string = A/B_C.dat -- sorts by B (alphabetical) then C (numerical)
// Done to sort negative numbers 
struct cmp_isov{
	bool operator()(const std::string& lhs, const std::string& rhs){
                if (lhs==rhs) return false;
                std::string comp = "/";
                std::vector<std::string> spllhs = Split(lhs,comp);
		std::vector<std::string> splrhs = Split(rhs,comp);
                int entryL = spllhs.size();
                int entryR = splrhs.size();
                std::string comp3 = "_";
                spllhs = Split(spllhs[entryL-1],comp3);
		splrhs = Split(splrhs[entryR-1],comp3);
                entryL = spllhs.size();
                entryR = splrhs.size();
		std::string comp2 = ".dat";
                std::vector<std::string> spllhsR = Split(spllhs[entryL-1],comp2);
                std::vector<std::string> splrhsR = Split(splrhs[entryR-1],comp2);
                if(spllhs[0]==splrhs[0])
		    return atof(spllhsR[0].c_str())<atof(splrhsR[0].c_str());
		else
		    return spllhs[0]<splrhs[0];
	}
};
//=============================================================================
#endif
//=============================================================================
