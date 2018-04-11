//=============================================================================
// Program for outputting the net yields for a set of masses and metallicities
//=============================================================================
#include "utils_ch.h"
#include "utils.h"
#include "sfr.h"
#include "imf.h"
#include "ages.h"
#include "iarates.h"
#include "params.h"
#include "yields.h"
#include "grid.h"
#include <iomanip>
//=============================================================================
#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP
//=============================================================================
struct result{
	std::string element;
	int number;
	double yield;
};
struct by_number{
	bool operator()(result const &a, result const &b){
		return a.number<b.number;
	}
};
//=============================================================================
int main(int argc, char const *argv[])
{
	ModelParameters M(argv[1]);
	YieldsSet Y(M);
	VecDoub MetalGrid = {1e-4,5e-4,1e-3,5e-3,1e-2,5e-2};
	VecDoub MassGrid = {1.,2.,5.,8.,10.,20.,50.,100.};
	std::cout<<"Z M Element Yield\n";
	for(auto z: MetalGrid)
		for(auto m: MassGrid){
			std::vector<result> results;
			for(auto e:element_index)
				results.push_back({e.first,e.second,Y.yield(e.first,m,z)});
			std::sort(results.begin(),results.end(),by_number());
			for(auto r:results)
				std::cout<< std::setprecision(5)<<z<<" "<<m<<" "<<r.element<<" "<<r.yield<<std::endl;
		}
	return 0;
}
//=============================================================================
