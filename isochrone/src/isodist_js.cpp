#include <Python.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>
#include <vector>

#include "utils.h"
#include "boost_python_interface.h"
#include "GSLInterface/GSLInterface.h"
#include "coordtransforms.h"

#include "utils_iso.h"
#include "extinction.h"
#include "iso_grid.h"
#include "distances.h"

// ============================================================================

namespace distance_compute{

    // Initialization and setup functions

    std::unique_ptr<isochrone_grid<isochrone_johnson>> iso;
    std::unique_ptr<isochrone_grid<isochrone_padova>> iso_Padova;
    std::unique_ptr<isochrone_grid<isochrone_dartmouth>> iso_Dartmouth;
    std::unique_ptr<DistanceCalculator<isochrone_johnson>> isoD;
    std::unique_ptr<DistanceCalculator<isochrone_padova>> isoD_Padova;
    std::unique_ptr<DistanceCalculator<isochrone_dartmouth>> isoD_Dartmouth;

    std::shared_ptr<extinction_law> EL;
    std::shared_ptr<extinction_map> EM;

    std::unique_ptr<galaxy_prior> prior;

    void init_isochrone(const std::string & which,int thin=1,double feh_err=1e-5,double thin_mag=-1.){
        if(which=="BaSTI"){
            iso = make_unique<isochrone_grid<isochrone_johnson>>(which,thin,feh_err,thin_mag);
            isoD = make_unique<DistanceCalculator<isochrone_johnson>>(&*iso);
        }
        else if(which=="Padova"){
            iso_Padova = make_unique<isochrone_grid<isochrone_padova>>(which,thin,feh_err,thin_mag);
            isoD_Padova = make_unique<DistanceCalculator<isochrone_padova>>(&*iso_Padova);
        }
        else if(which=="Dartmouth"){
            iso_Dartmouth = make_unique<isochrone_grid<isochrone_dartmouth>>(which,thin,feh_err,thin_mag);
            isoD_Dartmouth = make_unique<DistanceCalculator<isochrone_dartmouth>>(&*iso_Dartmouth);
        }
        else if(which=="All"){
            iso = make_unique<isochrone_grid<isochrone_johnson>>("BaSTI",thin,feh_err,thin_mag);
            iso_Padova = make_unique<isochrone_grid<isochrone_padova>>("Padova",thin,feh_err,thin_mag);
            iso_Dartmouth = make_unique<isochrone_grid<isochrone_dartmouth>>("Dartmouth",thin,feh_err,thin_mag);
            isoD = make_unique<DistanceCalculator<isochrone_johnson>>(&*iso);
            isoD_Padova = make_unique<DistanceCalculator<isochrone_padova>>(&*iso_Padova);
            isoD_Dartmouth = make_unique<DistanceCalculator<isochrone_dartmouth>>(&*iso_Dartmouth);
        }
        EL = std::make_shared<schlafly2017_extinction_law>();
    }

    void load_emap(void){
    	EM = std::make_shared<sfd_extinction_map>(&*EL,
                                                  conv::StandardSolarPAUL);
    }

    void load_prior(std::string prior_type="2018",
                    VecDoub SolarPosition=conv::StandardSolarPAUL){
        std::cout<<prior_type<<std::endl;
        if(prior_type=="Binney"){
            prior = make_unique<binney_prior>(SolarPosition);
            std::cout<<"Binney 2014 prior loaded for isochrone fitting\n";
        }
        else if(prior_type=="2018"){
            prior = make_unique<new_prior_2018>(SolarPosition);
            std::cout<<"New 2018 prior loaded for isochrone fitting\n";
        }
        else if(prior_type=="2018_broad"){
            prior = make_unique<new_prior_2018_broad_age>(SolarPosition);
            std::cout<<"New 2018 broad age prior loaded for isochrone fitting\n";
        }
        else if(prior_type=="2018_flat"){
            prior = make_unique<new_prior_2018_flat_age>(SolarPosition);
            std::cout<<"New 2018 flat age prior loaded for isochrone fitting\n";
        }
        else{
            prior = make_unique<galaxy_prior>(SolarPosition);
            std::cout<<"Flat prior loaded for isochrone fitting\n";
        }
    }

    double eval_prior(VecDoub X, double Z, double tau){
        return prior->prior(X,Z,tau);
    }

// ============================================================================
    // Sampling

    double photometric_distance_zero_age_pdf(double mu, VecDoub ZM,
                                             VecDoub mag, VecDoub inputs,
                                             VecDoub mag_errors,
                                  bool bprior,
                                  const std::string & which="BaSTI",
                                  std::vector<std::string> magstr={"H","J","K"}
){
        prior->set_prior(bprior);
        double l=inputs[0]; double b=inputs[1];
        double s = pow(10.,0.2*mu-2.);

        if(which=="BaSTI"){
            return isoD->photometric_distance_zero_age_pdf({l,b,s},ZM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else if(which=="Padova"){
            return isoD_Padova->photometric_distance_zero_age_pdf({l,b,s},ZM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else if(which=="Dartmouth"){
            return isoD_Dartmouth->photometric_distance_zero_age_pdf({l,b,s},ZM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else{
            std::cerr<<"Isochrone type not valid"<<std::endl;
            return 0.;
        }
    }
    double photometric_distance_pdf(double mu, VecDoub ZAM, VecDoub mag, VecDoub inputs,VecDoub mag_errors,
                                  bool bprior,
                                  const std::string & which="BaSTI",
                                  std::vector<std::string> magstr={"H","J","K"}
){
        prior->set_prior(bprior);
        double l=inputs[0]; double b=inputs[1];
        double s = pow(10.,0.2*mu-2.);

        if(which=="BaSTI"){
            return isoD->photometric_distance_pdf({l,b,s},ZAM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else if(which=="Padova"){
            return isoD_Padova->photometric_distance_pdf({l,b,s},ZAM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else if(which=="Dartmouth"){
            return isoD_Dartmouth->photometric_distance_pdf({l,b,s},ZAM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else{
            std::cerr<<"Isochrone type not valid"<<std::endl;
            return 0.;
        }
    }
    double distance_pdf(double mu, VecDoub ZTM, VecDoub mag, VecDoub inputs,
                                          VecDoub mag_errors, VecDoub icov,
                                  bool bprior,
                                  const std::string & which="BaSTI",
                                  std::vector<std::string> magstr={"H","J","K"}
){
        prior->set_prior(bprior);
        double Z=inputs[0]; double Teff=inputs[1];
        double logg=inputs[2]; double l=inputs[3]; double b=inputs[4];

        // Either pass 3D error vector 1/(ERR_Z,ERR_TEFF,ERR_LOGG) in icov
        // or 6D vector of C_00,C_10,C_20,C_11,C_21,C_22 where C is the inverse of the covariance matrix
        // of Z,TEFF,LOGG
	bool alpha = (inputs.size()==6);

        if(icov.size()==3+alpha){
            icov.push_back(icov[1]*icov[1]);
            icov.push_back(0.);
            if(alpha) icov.push_back(0.);
            icov.push_back(icov[2]*icov[2]);
            if(alpha){
        		icov.push_back(0.);
        		icov.push_back(icov[3]*icov[3]);
    	    }
            icov[0]*=icov[0];icov[1]=0.;icov[2]=0.;
            if(alpha) icov[3]=0.;
    	}

        double s = pow(10.,0.2*mu-2.);
	VecDoub ztg = {Z,Teff,logg};
        if(alpha) ztg.push_back(inputs[5]);

	if(which=="BaSTI"){
            return isoD->distance_pdf(ztg,{l,b,s},icov,ZTM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else if(which=="Padova"){
            return isoD_Padova->distance_pdf(ztg,{l,b,s},icov,ZTM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else if(which=="Dartmouth"){
            return isoD_Dartmouth->distance_pdf(ztg,{l,b,s},icov,ZTM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else{
            std::cerr<<"Isochrone type not valid"<<std::endl;
            return 0.;
        }
    }
    double distance_pdf_es(double mu, VecDoub ZTM, VecDoub mag, VecDoub inputs,
                                          VecDoub mag_errors, VecDoub icov,
                                  bool bprior,
                                  const std::string & which="BaSTI",
                                  std::vector<std::string> magstr={"H","J","K"}
){
        prior->set_prior(bprior);
        double Z=inputs[0]; double Teff=inputs[1];
        double logg=inputs[2]; double l=inputs[3]; double b=inputs[4];

        // Either pass 3D error vector 1/(ERR_Z,ERR_TEFF,ERR_LOGG) in icov
        // or 6D vector of C_00,C_10,C_20,C_11,C_21,C_22 where C is the inverse of the covariance matrix
        // of Z,TEFF,LOGG
        bool alpha = (inputs.size()==6);

        if(icov.size()==3+alpha){
            icov.push_back(icov[1]*icov[1]);
            icov.push_back(0.);
            if(alpha) icov.push_back(0.);
            icov.push_back(icov[2]*icov[2]);
            if(alpha){
                icov.push_back(0.);
                icov.push_back(icov[3]*icov[3]);
            }
            icov[0]*=icov[0];icov[1]=0.;icov[2]=0.;
            if(alpha) icov[3]=0.;
        }

        double s = pow(10.,0.2*mu-2.);
        VecDoub ztg = {Z,Teff,logg};
        if(alpha) ztg.push_back(inputs[5]);

        if(which=="BaSTI"){
            return isoD->distance_pdf_es(ztg,{l,b,s},icov,ZTM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else if(which=="Padova"){
            return isoD_Padova->distance_pdf_es(ztg,{l,b,s},icov,ZTM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else if(which=="Dartmouth"){
            return isoD_Dartmouth->distance_pdf_es(ztg,{l,b,s},icov,ZTM,mag,mag_errors,magstr,prior.get(),EM);
        }
        else{
            std::cerr<<"Isochrone type not valid"<<std::endl;
            return 0.;
        }
    }


// ============================================================================
    // Computing full distance pdf

    VecDoub prob_distance(VecDoub mag, VecDoub inputs,
                          VecDoub mag_errors, VecDoub covar,
                          bool bprior,
                          const std::string & which="BaSTI",
                          std::vector<std::string> magstr={"H","J","K"},
                          bool extinction=false, double Aprior=0.,
                          bool return_pdf=false,
                          double parallax=0.,
                          double parallax_error=-1.,
                          double mass=0.,
                          double mass_error=-1.){
        prior->set_prior(bprior);
        double Z=inputs[0]; double Teff=inputs[1];
        double logg=inputs[2]; double l=inputs[3]; double b=inputs[4];

	// Either pass 3D error vector (ERR_Z,ERR_TEFF,ERR_LOGG) in covar
	// or 6D vector of C_00,C_10,C_20,C_11,C_21,C_22 where C is the covariance matrix
	// of Z,TEFF,LOGG

	if(covar.size()==3){
		covar.push_back(covar[1]*covar[1]);
		covar.push_back(0.);
		covar.push_back(covar[2]*covar[2]);
		covar[0]*=covar[0];covar[1]=0.;covar[2]=0.;
	}

        double err_Z=sqrt(covar[0]); double err_Teff=sqrt(covar[3]);
        double err_logg=sqrt(covar[5]);

        if(which=="BaSTI"){
            if(extinction) return isoD->prob_distance_extinct(mag,Z,Teff,logg,l,b,mag_errors,covar,prior.get(),return_pdf,5.,magstr,&*EM,Aprior,parallax,parallax_error);
            else return isoD->prob_distance(mag,Z,Teff,logg,l,b,
                                            mag_errors,covar,prior.get(),
                                            return_pdf,5.,magstr,
                                            parallax,parallax_error,
                                            mass,mass_error);
        }
        else if(which=="Padova"){
            if(extinction) return isoD_Padova->prob_distance_extinct(mag,Z,Teff,logg,l,b,mag_errors,covar,prior.get(),return_pdf,5.,magstr,&*EM,Aprior,parallax,parallax_error);
            else return isoD_Padova->prob_distance(mag,Z,Teff,logg,l,b,
                                                   mag_errors,
                                                   covar,prior.get(),
                                                   return_pdf,5.,magstr,
                                                   parallax,parallax_error,
                                                   mass,mass_error);
        }
        else if(which=="Dartmouth"){
            if(extinction) return isoD_Dartmouth->prob_distance_extinct(mag,Z,Teff,logg,l,b,mag_errors,covar,prior.get(),return_pdf,5.,magstr,&*EM,Aprior,parallax,parallax_error);
            else return isoD_Dartmouth->prob_distance(mag,Z,Teff,logg,l,b,
                                                      mag_errors,
                                                      covar,prior.get(),
                                                      return_pdf,5.,magstr,
                                                      parallax,parallax_error,
                                                      mass,mass_error);
        }
        else{
            std::cerr<<"Isochrone type not valid"<<std::endl;
            return {0.,0.,0.};
        }
    }

    VecDoub prob_distance_alpha(VecDoub mag, VecDoub inputs,
                                VecDoub mag_errors, VecDoub covar,
                                bool bprior,
                                const std::string & which="BaSTI",
                                std::vector<std::string> magstr={"H","J","K"},
                                bool extinction=false,
                                double Aprior=0., bool return_pdf=false){
        prior->set_prior(bprior);
        double Z=inputs[0]; double Teff=inputs[1];
        double logg=inputs[2]; double alpha=inputs[3]; double l=inputs[4]; double b=inputs[5];

    // Either pass 4D error vector (ERR_Z,ERR_TEFF,ERR_LOGG,ERR_ALP) in covar
    // or 10D vector of C_00,C_10,C_20,C_30,C_11,C_21,C_31,C_22,C_32,C_33 where C is the covariance matrix
    // of Z,TEFF,LOGG

    if(covar.size()==4){
        covar.push_back(covar[1]*covar[1]);
        covar.push_back(0.);
        covar.push_back(0.);
        covar.push_back(covar[2]*covar[2]);
        covar.push_back(0.);
        covar.push_back(covar[3]*covar[3]);
        covar[0]*=covar[0];covar[1]=0.;covar[2]=0.;covar[3]=0.;
    }

        double err_Z=sqrt(covar[0]); double err_Teff=sqrt(covar[4]);
        double err_logg=sqrt(covar[7]);
        double err_alpha=sqrt(covar[9]);

        if(which=="BaSTI"){
            if(extinction) return isoD->prob_distance_extinct(mag,Z,Teff,logg,l,b,mag_errors,covar,prior.get(),return_pdf,5.,magstr,&*EM,Aprior);
            else return isoD->prob_distance_alpha(mag,Z,Teff,logg,alpha,l,b,mag_errors,covar,prior.get(),5.,magstr);
        }
        else if(which=="Padova"){
            if(extinction) return isoD_Padova->prob_distance_extinct(mag,Z,Teff,logg,l,b,mag_errors,covar,prior.get(),return_pdf,5.,magstr,&*EM,Aprior);
            else return isoD_Padova->prob_distance_alpha(mag,Z,Teff,logg,alpha,l,b,mag_errors,covar,prior.get(),5.,magstr);
        }
        else if(which=="Dartmouth"){
            if(extinction) return isoD_Dartmouth->prob_distance_extinct(mag,Z,Teff,logg,l,b,mag_errors,covar,prior.get(),return_pdf,5.,magstr,&*EM,Aprior);
            else return isoD_Dartmouth->prob_distance_alpha(mag,Z,Teff,logg,alpha,l,b,mag_errors,covar,prior.get(),5.,magstr);
        }
        else{
            std::cerr<<"Isochrone type not valid"<<std::endl;
            return {0.,0.,0.};
        }
    }
    VecDoub prob_distance_extinctprior(VecDoub mag, VecDoub inputs,
                                       VecDoub mag_errors, VecDoub covar,
                                       bool bprior,VecDoub Aprior,
                                       VecDoub sigAprior, VecDoub log_dist_AV,
                                       const std::string & which="BaSTI",
                                  std::vector<std::string> magstr={"H","J","K"},
                                  bool return_pdf=false,
                                  double parallax=0.,
                                  double parallax_error=-1.,
                                  double mass=0.,
                                  double mass_error=-1.){
        prior->set_prior(bprior);
        double Z=inputs[0]; double Teff=inputs[1];
        double logg=inputs[2]; double l=inputs[3]; double b=inputs[4];

    // Here we pass in a prior on log(A_V), the associated spread and the log(distances) at which they are evaluated.

    // Either pass 3D error vector (ERR_Z,ERR_TEFF,ERR_LOGG) in covar
    // or 6D vector of C_00,C_10,C_20,C_11,C_21,C_22 where C is the covariance matrix
    // of Z,TEFF,LOGG

    if(covar.size()==3){
        covar.push_back(covar[1]*covar[1]);
        covar.push_back(0.);
        covar.push_back(covar[2]*covar[2]);
        covar[0]*=covar[0];covar[1]=0.;covar[2]=0.;
    }

        double err_Z=sqrt(covar[0]); double err_Teff=sqrt(covar[3]);
        double err_logg=sqrt(covar[5]);

        if(which=="BaSTI"){
            return isoD->prob_distance_extinctprior(mag,Z,Teff,logg,l,b,mag_errors,covar,prior.get(),return_pdf,5.,magstr,Aprior,sigAprior,log_dist_AV,&*EL,parallax,parallax_error,
                                                      mass,mass_error);
        }
        else if(which=="Padova"){
            return isoD_Padova->prob_distance_extinctprior(mag,Z,Teff,logg,l,b,mag_errors,covar,prior.get(),return_pdf,5.,magstr,Aprior,sigAprior,log_dist_AV,EL.get(),parallax,parallax_error,
                                                      mass,mass_error);
        }
        else if(which=="Dartmouth"){
            return isoD_Dartmouth->prob_distance_extinctprior(mag,Z,Teff,logg,l,b,mag_errors,covar,prior.get(),return_pdf,5.,magstr,Aprior,sigAprior,log_dist_AV,&*EL,parallax,parallax_error,
                                                      mass,mass_error);
        }
        else{
            std::cerr<<"Isochrone type not valid"<<std::endl;
            return {0.,0.,0.};
        }
    }

// ============================================================================
    // Auxiliary functions

    double get_es(double Z, double age, double M, std::string which){
        if(which=="BaSTI"){
            return age/iso->max_age(Z,M);
        }
        else if(which=="Padova"){
            return age/iso_Padova->max_age(Z,M);
        }
        else if(which=="Dartmouth"){
            return age/iso_Dartmouth->max_age(Z,M);
        }
        else{
            std::cerr<<"Invalid isochrone"<<std::endl;
            return -9999.;
        }
    }
    double get_magnitude(double Z, double age, double M, std::string mag, std::string which,bool interp){
        if(which=="BaSTI"){
            if(interp) return iso->interp_es(Z,age,M,{mag})[2];
            std::vector<int> near = iso->find_nearest(Z,age,M);
    	    return iso->iso(near[0],near[1])->mag(near[2],mag);
    	}
        else if(which=="Padova"){
            if(interp) return iso_Padova->interp_es(Z,age,M,{mag})[2];
            std::vector<int> near = iso_Padova->find_nearest(Z,age,M);
    	    return iso_Padova->iso(near[0],near[1])->mag(near[2],mag);
    	}
        else if(which=="Dartmouth"){
            if(interp) return iso_Dartmouth->interp_es(Z,age,M,{mag})[2];
            std::vector<int> near = iso_Dartmouth->find_nearest(Z,age,M);
    	    return iso_Dartmouth->iso(near[0],near[1])->mag(near[2],mag);
    	}
    	else{
    		std::cerr<<"Invalid isochrone"<<std::endl;
    		return -9999.;
    	}
    }
    VecDoub get_teff_logg(double Z, double age, double M, std::string which,bool interp, double alpha=0.){
        double metal = Z+log10(0.638*pow(10.,alpha)+0.362);
        if(which=="BaSTI"){
            if(interp){
                auto xx = iso->interp_es(metal,age,M,{"J"});
                return {xx[0],xx[1]};
            }
            std::vector<int> near = iso->find_nearest(metal,age,M);
            return {iso->iso(near[0],near[1])->logTeff(near[2]),
                    iso->iso(near[0],near[1])->logg(near[2])};
        }
        else if(which=="Padova"){
            if(interp){
                auto xx = iso_Padova->interp_es(metal,age,M,{"J"});
                return {xx[0],xx[1]};
            }
            std::vector<int> near = iso_Padova->find_nearest(metal,age,M);
            return {iso_Padova->iso(near[0],near[1])->logTeff(near[2]),
                    iso_Padova->iso(near[0],near[1])->logg(near[2])};
        }
        else if(which=="Dartmouth"){
            if(interp){
                auto xx = iso_Dartmouth->interp_es(metal,age,M,{"J"});
                return {xx[0],xx[1]};
            }
            std::vector<int> near = iso_Dartmouth->find_nearest(metal,age,M);
            return {iso_Dartmouth->iso(near[0],near[1])->logTeff(near[2]),
                    iso_Dartmouth->iso(near[0],near[1])->logg(near[2])};
        }
        else{
            std::cerr<<"Invalid isochrone"<<std::endl;
            return {-9999.,-9999.};
        }
    }
}

// ============================================================================

BOOST_PYTHON_MODULE(isodist_js)
{

    docstring_options doc_options(true);
    numeric::array::set_module_and_type("numpy", "ndarray");

    // Initialization functions
    def("init_isochrone", distance_compute::init_isochrone);
    def("load_emap", distance_compute::load_emap);
    def("load_prior", distance_compute::load_prior);

    // For evaluating pdfs
    def("distance_pdf", distance_compute::distance_pdf);
    def("distance_pdf_es", distance_compute::distance_pdf_es);
    def("photometric_distance_zero_age_pdf", distance_compute::photometric_distance_zero_age_pdf);
    def("photometric_distance_pdf", distance_compute::photometric_distance_pdf);

    // For constructing 1d distance pdf
    def("prob_distance", distance_compute::prob_distance);
    def("prob_distance_alpha", distance_compute::prob_distance_alpha);
    def("prob_distance_extinctprior",
        distance_compute::prob_distance_extinctprior);

    // Auxiliary functions
    def("get_es", distance_compute::get_es);
    def("get_magnitude", distance_compute::get_magnitude);
    def("get_teff_logg", distance_compute::get_teff_logg);

    // Prior
    def("prior", distance_compute::eval_prior);

    to_python_converter<VecDoub, vector_to_ndarray<double>>();
    vector_from_ndarray<double>();
    to_python_converter<std::vector<std::string>, vector_to_ndarray<std::string>>();
    vector_from_ndarray<std::string>();
    import_array();

}
