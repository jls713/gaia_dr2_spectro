#include "params.h"
using namespace H5;
//=============================================================================
// Constructor
ModelParameters::ModelParameters(std::string filename){
	std::ifstream inFile(filename);
	if(!inFile.good())
		LOG(ERROR)<<"Cannot open "<<filename<<std::endl;
	inFile>>parameters;
	inFile.close();
}
//=============================================================================
// Printing and output
void ModelParameters::print(void){
	std::cout<<parameters<<std::endl;
}
void ModelParameters::pretty_print(std::ostream& out){

	out<<"\nModel Parameters\n==============\n";

	out<<"Fundamentals\n-------------\n";
	out<<"Galaxy Age="
	   <<parameters["fundamentals"]["GalaxyAge"]<<std::endl;
	out<<"Minimum mass ="
	   <<parameters["fundamentals"]["MinimumMass"]<<std::endl;
	out<<"Maximum mass ="
	   <<parameters["fundamentals"]["MaximumMass"]<<std::endl;
	out<<"IMF="
	   <<parameters["fundamentals"]["IMF"]<<std::endl;
	out<<"SFR="
	   <<parameters["fundamentals"]["SFR"]<<std::endl;
	out<<"Stellar ages="
	   <<parameters["fundamentals"]["lifetimes"]<<std::endl;
	out<<"-------------\n";

	out<<"Grids\n-------------\n";
	out<<"Minimum radius ="
	   <<parameters["grids"]["MinimumRadius"]<<std::endl;
	out<<"Maximum radius ="
	   <<parameters["grids"]["MaximumRadius"]<<std::endl;
	out<<"Number of radial grid points="
	   <<parameters["grids"]["RadialGridPoints"]<<std::endl;
	out<<"Number of time grid points="
	   <<parameters["grids"]["AgeGridPoints"]<<std::endl;
	out<<"-------------\n";

	out<<"Yields\n-------------\n";
	out<<"AGB yields="
	   <<parameters["yields"]["AGB"]<<std::endl;
	out<<"Super AGB yields="
	   <<parameters["yields"]["SuperAGB"]<<std::endl;
	out<<"Type II yields="
	   <<parameters["yields"]["typeII"]<<std::endl;
	out<<"Type Ia yields="
	   <<parameters["yields"]["typeIa"]<<std::endl;
	out<<"-------------\n";

	out<<"Type Ia Rates\n-------------\n";
	out<<"Form ="
	   <<parameters["typeIa"]["Form"]<<std::endl;
	out<<"Minimum Binary Mass ="
	   <<parameters["typeIa"]["MinimumBinaryMass"]<<std::endl;
	out<<"Maximum Binary Mass ="
	   <<parameters["typeIa"]["MaximumBinaryMass"]<<std::endl;
	out<<"Binary Fraction ="
	   <<parameters["typeIa"]["BinaryFraction"]<<std::endl;
	out<<"-------------\n";

	out<<"Elements = {";
	auto elements_list = parameters["elements"];
	for(auto i:elements_list)
		std::cout<<i<<",";
	out<<"}\n";
	out<<"-------------\n";

	out<<std::endl;
}
void ModelParameters::write_hdf5(H5File &fout){
	auto w = parameters.dump().c_str();
	auto s = parameters.dump();
	StrType dtype(PredType::C_S1,s.length()+1);
	DataSpace dspace(H5S_SCALAR);
	DataSet dataset = fout.createDataSet("parameters",dtype,dspace);
	dataset.write(w,dtype);
	return;
}
//=============================================================================
