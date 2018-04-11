#include "mass_function.h"

MassFunction::MassFunction(json parameters){
	IMFflag = parameters["MassFunction"]["IMFflag"];
	ModelParameters mparameters(parameters);
	mparameters.parameters["fundamentals"]["MinimumMass"] = parameters["MassFunction"]["MinMass"];
	mparameters.parameters["fundamentals"]["MaximumMass"] = parameters["MassFunction"]["MaxMass"];
	mparameters.parameters["data_folder"] = parameters["MassFunction"]["LifetimeFolder"];
	IMF = imf_types[parameters["MassFunction"]["IMF"]](mparameters);
	Lifetime = life_types[parameters["MassFunction"]["Lifetime"]](mparameters);
	double MinMass = parameters["MassFunction"]["MinMass"];
	double MaxMass = parameters["MassFunction"]["MaxMass"];
	LogMassGrid = create_range(log(MinMass),
	                           log(MaxMass),
	                           NM);
	double deltaLogM = LogMassGrid[1]-LogMassGrid[0];
	CMFgrid.push_back(0.);
	for(unsigned m=1;m<LogMassGrid.size();++m){
		double mass = exp(LogMassGrid[m]-deltaLogM);
		CMFgrid.push_back(deltaLogM*mass*(*IMF)(mass)+CMFgrid[m-1]);
	}
}

double MassFunction::CMF(double MaxMass){
	return linterp(LogMassGrid,CMFgrid,log(MaxMass),"const");
}
