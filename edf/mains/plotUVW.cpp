#include "edf.h"

int main(int argc, char *argv[]){

	std::string galpot_file=parameters()["potential"];
	GalPot Pot(galpot_file);
	std::string action_file=parameters()["actions"];
	Actions_AxisymmetricFudge_InterpTables Tab(&Pot,action_file);

	sb15_edf edf(&Pot,&Tab);

	std::string params_file=parameters()["edf_params"];
	edf.readParams(params_file);

	std::ofstream outFile; outFile.open(argv[1]);
	double IE=1e-2; double RSolar = conv::StandardSolar[0];
	for(double V=-150.;V<150.;V+=5.){
	        outFile<<V<<" "<<V-170.<<" "<<V<<" "<<edf.UHist(RSolar,0.,V,IE)<<" "<<edf.VHist(RSolar,0.,-V+170.,IE)<<" "<<edf.WHist(RSolar,0.,V,IE)<<std::endl;
	}
	outFile.close();
}
