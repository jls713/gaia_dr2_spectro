#include "edf.h"

int main(int argc, char *argv[]){

	std::string galpot_file=edf_parameters()["potential"];
	std::cout<<galpot_file<<std::endl;
	GalPot Pot(galpot_file);
	std::string action_file=edf_parameters()["actions"]["dir"];
	Actions_AxisymmetricFudge_InterpTables Tab(&Pot,action_file);

	VecDoub solar_motion = edf_parameters()["solar_motion"];

	sb15_edf edf(&Pot,&Tab,{0.,0.},solar_motion);

	std::string params_file=edf_parameters()["edf_params"];
	edf.readParams(params_file);

	auto F = edf_parameters();
	if (F.find("edf_json_params") != F.end()){
		std::string flname=F["edf_json_params"];
		edf.setParams(flname,true);
	}

	std::ofstream outFile; outFile.open(argv[1]);
	double IE=1e-2; double RSolar = solar_motion[0], zsolar = solar_motion[1];
	for(double V=-1.;V<.6;V+=0.01){
	        outFile<<V<<" "<<edf.FHist(RSolar,zsolar,V,IE)<<std::endl;
	}
	outFile.close();
}
