#include "edf.h"

int main(){
	double FSplit = -0.2;

	std::string galpot_file=edf_parameters()["potential"];
	GalPot Pot(galpot_file);
	std::string action_file=edf_parameters()["potential"];
	Actions_AxisymmetricFudge_InterpTables Tab(&Pot,action_file);

	sb15_edf edf(&Pot,&Tab);

	std::string params_file=edf_parameters()["edf_params"];
	edf.readParams(params_file);

	edf.TurnOffHalo();
	edf.printParams();

	VecDoub StandardSolar = edf.SunCoords();

	double IE = 1e-2;
	std::ofstream outFile; outFile.open("global_profiles/densityR_metalsplit.dat");
	for(double R=0.1;R<16.;R+=0.4){
		edf.TurnOnThin();edf.TurnOffThick();
		outFile<<R<<" "<<edf.density(R,0.,IE,{-10.,FSplit})<<" "<<edf.density(R,0.,IE,{FSplit,10.})<<" ";
		edf.TurnOnThick();edf.TurnOffThin();
		outFile<<edf.density(R,0.,IE,{-10.,FSplit})<<" "<<edf.density(R,0.,IE,{FSplit,10.})<<std::endl;
	}
	outFile.close();

	std::ofstream outFile2; outFile2.open("global_profiles/densityz_metalsplit.dat");
	for(double z=0.;z<4.6;z+=0.1){
		edf.TurnOnThin();edf.TurnOffThick();
		outFile2<<z<<" "<<edf.density(StandardSolar[0],z,IE,{-10.,FSplit})<<" "<<edf.density(StandardSolar[0],z,IE,{FSplit,10.})<<" ";
		edf.TurnOnThick();edf.TurnOffThin();
		outFile2<<edf.density(StandardSolar[0],z,IE,{-10.,FSplit})<<" "<<edf.density(StandardSolar[0],z,IE,{FSplit,10.})<<std::endl;
	}
	outFile2.close();
}
