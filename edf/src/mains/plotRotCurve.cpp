
#include "utils.h"
#include "potential.h"

int main(int argc, char const *argv[])
{
	if(argc<2){
		std::cerr<<"Must pass Tpot file as second parameter\n";
		return 0;
	}
	GalPot G(argv[1]);
	for(auto i: create_range(.1,30.,200))
		std::cout<<i<<" "<<G.Vc(i)<<std::endl;
	std::cerr<<G.Vc(8.)<<std::endl;
	return 0;
}
