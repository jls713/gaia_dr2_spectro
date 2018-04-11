#ifndef GR_
#define GR_

#include "utils.h"
#include "edf.h"
#include <cmath>
#include <iostream>

class GRdata{
	public:
		VecDoub GR_ZZ, GR_GR, GR_Err1, GR_Err2;
	GRdata(){
		std::ifstream inFileGR; inFileGR.open("data/Gilmore_Reid/GR.dat");
		if(!inFileGR.is_open()){std::cerr<<"GR file won't open."<<std::endl;}
		else std::cerr<<"GR file loaded."<<std::endl;
		double ZZ, GR, Err1, Err2;
		while(inFileGR>>ZZ>>GR>>Err1>>Err2){
			GR_ZZ.push_back(ZZ);GR_GR.push_back(GR);
			GR_Err1.push_back(Err1);GR_Err2.push_back(Err2);
		}
		inFileGR.close();
	}
};

#endif
