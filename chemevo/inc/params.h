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
template<typename T>
std::string toString(const T& t) {
    return std::to_string(t);
}
std::string toString(const char* t);
std::string toString(const std::string& t);
//=============================================================================
template<typename T, typename T2>
std::vector<T2> extract_params(T F, std::vector<std::string> vars,
                               std::vector<T2> defaults, bool log=true){
    std::vector<T2> rslt;
    if(defaults.size()!=vars.size()){
        std::string varst="";
        std::string dst="";
        for(auto d:defaults) dst+=toString(d);
        for(auto v:vars) varst+=v;
        LOG(INFO)<<varst<<" "<<dst<<" different length";
        throw std::invalid_argument("Mismatching lengths "+varst+" "+dst);
    }
    for(unsigned ii=0;ii<vars.size();++ii){
        if (F.find(vars[ii]) == F.end()) {
            if(log)
            	LOG(INFO)<<"'"<<vars[ii]<<"'"
	            <<" not found in parameters file. Setting default for '"
	            <<vars[ii]<<"' = "<<defaults[ii];
	        rslt.push_back(defaults[ii]);
        }
        else
	        rslt.push_back(F[vars[ii]]);
    }
    return rslt;
}
template<typename T, typename T2>
T2 extract_param(T F, std::string var,
                  T2 defaultP, bool log=true){
	std::vector<std::string> vars={var};
	std::vector<T2> defaults={defaultP};
    return extract_params(F, vars, defaults, log)[0];
}
template<typename T>
int check_param_given(T F, std::string entry, bool log=true){
	if (F.find(entry) == F.end()) {
		if (log)
			LOG(INFO)<<"No "<<entry<<" parameter found in parameters file\n";
		return 1;
	}
	return 0;
}
template<typename T1, typename T2>
int check_param_matches_list(T1 F, T2 types_dict, std::string varname){
	if(types_dict.count(F)!=1){
		std::string types = "";
		for(auto i:types_dict) types+=i.first+",";
		LOG(INFO)<<"Invalid "<<varname<<": "<<F<<std::endl;
		LOG(INFO)<<"Options are:"<<types<<std::endl;
		return 1;
	}
	return 0;
}
template<typename T1, typename T2>
int check_param_given_matches_list(T1 F, T2 types_dict, std::string varname){
	if(check_param_given(F, varname)){
		return 1;
	}
	else return check_param_matches_list(F[varname], types_dict, varname);
}
template<typename T1, typename T2>
int check_param_given_matches_list_form(T1 F, T2 types_dict,
                                        std::string varname){
	if(check_param_given(F, varname)){
		return 1;
	}
	else return check_param_matches_list(F[varname]["Form"],
	                                     types_dict, varname);
}
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
