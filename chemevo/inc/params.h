#ifndef PARAMS_H
#define PARAMS_H
//=============================================================================
#include <fstream>
//=============================================================================
#include "easylogging++.h"
#include "json.hpp"
#include "H5Cpp.h"
using json = nlohmann::json;
//=============================================================================
#include "utils_ch.h"
//=============================================================================
/**
 * @brief Model Parameters storage class
 * @details Uses functionality of json library to read a json parameter file
 * and store all the model parameters
 */
class ModelParameters{
private:
public:
	json parameters; // reads in parameters from json file
	// Constructor -- reads json file
	ModelParameters(std::string filename);
	ModelParameters(ModelParameters &p){parameters=p.parameters;}
	ModelParameters(json &p){parameters=p;}
	//=========================================================================
	// Prints and outputs
	void print(void);
	void pretty_print(std::ostream &out=std::cout);
	void output(std::string filename);
	void write_hdf5(H5::H5File &fout);
	//=========================================================================
	// Getters
	inline double min_typeII_SN_mass(void){
		return parameters["typeII"]["Min_typeII_SN_mass"];}
	inline double max_mass(void){
		return parameters["fundamentals"]["MaximumMass"];}
	inline double min_mass(void){
		return parameters["fundamentals"]["MinimumMass"];}
};
//=============================================================================
/**
 * @brief Unique pointer creator
 * @details Create a unique pointer of base class T and derived class C using arguments Args
 *
 * @param args arguments of constructor
 * @tparam T Base class
 * @tparam C Derived Class
 * @tparam typename ...Args argument types
 * @return pointer to base class
 */
template<typename T, typename C, typename ...Args>
std::unique_ptr<T> createInstance(Args ...args) {
	return make_unique<C>(args...);
}
/**
 * @brief Shared pointer creator
 * @details Create a shared pointer of base class T and derived class C using arguments Args
 *
 * @param args arguments of constructor
 * @tparam T Base class
 * @tparam C Derived Class
 * @tparam typename ...Args argument types
 * @return pointer to base class
 */
template<typename T,typename C, typename ...Args>
std::shared_ptr<T> createSharedInstance(Args ...args) {
	return std::make_shared<C>(args...);
}
// Shared map type -- map of strings and shared pointers
template<typename T, typename ...Args>
using shared_map
	= std::map<std::string, std::shared_ptr<T>(*)(Args ...args)>;
// Unique map type -- map of strings and unique pointers
template<typename T, typename ...Args>
using unique_map
	= std::map<std::string, std::unique_ptr<T>(*)(Args ...args)>;
#endif
//=============================================================================
