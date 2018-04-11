//=============================================================================
// Play around with a chemical evolution model
//=============================================================================
#include "model.h"
#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP
//=============================================================================
int main(int argc, char const *argv[])
{
    START_EASYLOGGINGPP(argc,argv);

    try{
	    Model M("example_params.json");
	    std::cout<<M.AGBEnrichmentRate("Fe",20.,4.94118)<<std::endl;
	}
	catch(std::exception const &e){
		LOG(INFO)<<e.what()<<std::endl;
		std::cerr<<e.what()<<std::endl;
	}
}
//=============================================================================
