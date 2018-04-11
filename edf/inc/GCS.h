#ifndef GCS_
#define GCS_

#include "utils.h"
#include "edf.h"
#include <cmath>
#include <iostream>

VecDoub find_max_values(sb15_edf_sf *edf, bool atSun);
VecDoub Sampler(double l, double b, sb15_edf_sf *edf, VecDoub Err, VecDoub max, bool atSun, std::ofstream &outFile);

class GCS_data{
	// For storing GCS data
	public:
		MatDoub InputCoords; // r, g-r, F, l, b, v_||, mu
		MatDoub InputCoords_e;
		MatDoub GalIn;
		VecDoub Metal, Dist;
		double tmp;
		MatDoub ActionTable;
		int NDATA, NPOINTS, ERROR_SAMPLES;

	GCS_data(std::string IN);
	void find_actions(Potential_JS *Pot, Actions_AxisymmetricFudge_InterpTables *Tab, int readyForLL);
	void make_fake_table(sb15_edf_sf *edf, bool atSun, std::string outfile);
	double LogLike(sb15_edf_sf *edf);
};

#endif
