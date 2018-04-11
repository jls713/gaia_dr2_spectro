#include "edf.h"
#include "GR.h"

int main(){

	std::string galpot_file=edf_parameters()["potential"];
	std::cout<<galpot_file<<std::endl;
	GalPot Pot(galpot_file);
	std::string action_file=edf_parameters()["actions"];
	Actions_AxisymmetricFudge_InterpTables Tab(&Pot,action_file);

	sb15_edf edf(&Pot,&Tab);

	std::string params_file=edf_parameters()["edf_params"];
	edf.readParams(params_file);

	edf.printParams();
	edf.TurnOnHalo();
	GRdata *Gil;
	Gil = new GRdata();
	double IE=1e-3; double RSolar = conv::StandardSolar[0];
	double Density=edf.density(RSolar,0., IE);
	std::cout<<"Density: "<<Density<<std::endl;
	std::ofstream outFile; outFile.open("data/Gilmore_Reid/modelGR.test");
	for(int i=0;i<Gil->GR_ZZ.size();i++){
	    outFile<<Gil->GR_ZZ[i]<<" "<<edf.density(RSolar,Gil->GR_ZZ[i],IE)/Density<<std::endl;
	}
	outFile.close();
}
