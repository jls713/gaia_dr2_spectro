#include "edf.h"

int main(int argc, char *argv[]){


	std::string galpot_file=edf_parameters()["potential"];
	GalPot Pot(galpot_file);
	std::string action_file=edf_parameters()["potential"];
	Actions_AxisymmetricFudge_InterpTables Tab(&Pot,action_file);

	sb15_edf edf(&Pot,&Tab);

	std::string params_file=edf_parameters()["edf_params"];
	edf.readParams(params_file);

	edf.printParams(std::cerr);

	edf.TurnOffThin();
	edf.TurnOffThick();
	edf.TurnOffHalo();

	double IE = 1e-3;
	std::ofstream outFile; outFile.open(argv[1]);
	for(double R=4.;R<=12.;R+=1.){
		edf.TurnOnThin();
		outFile<<R<<" "<<edf.density(R,1e-3,IE)<<" ";
		edf.TurnOffThin();
		edf.TurnOnThick();
		outFile<<edf.density(R,1e-3,IE)<<" ";
		edf.TurnOffThick();
		edf.TurnOnHalo();
		outFile<<edf.density(R,1e-3,IE)<<std::endl;
		edf.TurnOffHalo();
	}
	outFile.close();

	double R = conv::StandardSolar[0];
	if(argc>5) R = atof(argv[3]);

	std::ofstream outFile2; outFile2.open(argv[2]);
	for(double z=0.;z<4.6;z+=0.05){
		edf.TurnOnThin();
		outFile2<<z<<" "<<edf.density(R,z,IE)<<" ";
		edf.TurnOffThin();
		edf.TurnOnThick();
		outFile2<<edf.density(R,z,IE)<<" ";
		edf.TurnOffThick();
		edf.TurnOnHalo();
		outFile2<<edf.density(R,z,IE)<<std::endl;
		edf.TurnOffHalo();
	}
	outFile2.close();
}
