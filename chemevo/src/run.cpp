//=============================================================================
// Run a chemical evolution model
//=============================================================================
#include "model.h"
#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP
//=============================================================================
int main(int argc, char const *argv[])
{
    START_EASYLOGGINGPP(argc,argv);
	gsl_set_error_handler_off();
	setenv("CUBACORES","0",1);

    try{
	    std::string paramfile = "params/run_params.json";
	    if(argc>2) paramfile=std::string(argv[2]);
	    Model M(paramfile);
	    M.run();
	    std::string outputfile = "tmp.h5";
	    if(argc>1) outputfile=std::string(argv[1]);
	    M.write(outputfile);
	}
	catch(std::exception const &e){
	    LOG(INFO)<<e.what()<<std::endl;
	    std::cerr<<e.what()<<std::endl;
	}
}
//=============================================================================
