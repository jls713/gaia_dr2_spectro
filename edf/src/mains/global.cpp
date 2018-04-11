#include "edf.h"

int main(){

	std::string galpot_file=edf_parameters()["potential"];
	GalPot Pot(galpot_file);
	std::string action_file=edf_parameters()["potential"];
	Actions_AxisymmetricFudge_InterpTables Tab(&Pot,action_file);

	sb15_edf edf(&Pot,&Tab);

	std::string params_file=edf_parameters()["edf_params"];
	edf.readParams(params_file);

	std::ofstream outFile; outFile.open("GCSdata/globalRz.dat");
	double IE=1e-2;
	for(double R=2.;R<20.;R+=1.){
		for(double z=0.;z<.8;z+=.2){
	        outFile<<R<<" "<<z<<" "<<edf.density(R,z,IE)<<std::endl;
		}
	}
	outFile.close();
}
