#ifndef UTILS_CH
#define UTILS_CH
//=============================================================================
#include <memory>
#include <stdexcept>
#include <exception>
#include <string>
//=============================================================================
#include "utils.h"
#include "cuba.h"
#include "hdf5_reader.h"
#include "utils_iso.h"
//=============================================================================
const int nproc = 1;
const int SEED = time(0);
/**
 * @brief integrate using cuba
 * @details integrate general function <integrand> using cuba algorithm of type
 *
 * @param integrand integrand function
 * @param P additional parameters class
 * @param IE relative error required
 * @param AE absolute error required
 * @param type integration algorithm name ("Cuhre","Suave","Vegas","Divonne")
 * @param err returns error
 * @param str name of caller for error reporting
 * @return value of integral
 */
template<class c>
double integrate(integrand_t integrand, c *P, double IE, double AE, std::string type, double *err,std::string str){

    int neval,fail,nregions;
    double integral[1],error[1],prob[1];
    int NSIZE = P->x2min.size();
    double prod = 1.;
    for(int i=0;i<NSIZE;i++)prod*=(P->x2max[i]-P->x2min[i]);

    if(type=="Vegas")
        Vegas(NSIZE,nproc,integrand,P,1,IE,AE,0,SEED,
        MINEVAL,MAXEVAL,NSTART,NINCREASE,NBATCH,GRIDNO,STATEFILE,SPIN,
        &neval,&fail,integral,error,prob);

    else if (type=="Suave")
        Suave(NSIZE,nproc,integrand,P,1,IE,AE,0,SEED,
        MINEVAL,MAXEVAL,NNEW,FLATNESS,STATEFILE,SPIN,&nregions,
        &neval,&fail,integral,error,prob);

    else if (type=="Cuhre")
        Cuhre(NSIZE,nproc,integrand,P,1,IE,AE,0,
        MINEVAL, MAXEVAL, 0, STATEFILE,SPIN,
        &nregions, &neval, &fail, integral, error, prob);

    else
        Divonne(NSIZE,nproc,integrand,P,1,IE,AE,0,SEED,
        MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
        BORDER, MAXCHISQ, MINDEVIATION,
        NGIVEN, NSIZE, nullptr, NEXTRA, nullptr,STATEFILE,SPIN,
        &nregions, &neval, &fail, integral, error, prob);
    if(err)*err=prod*error[0];
    if(fail!=0)
        throw std::runtime_error("Error: Required accuracy not reached for "+str+".");
    return prod*integral[0];
}
//=============================================================================
template<class c>
void printProgBar( c percent ){
  std::string bar;

  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%     " << std::flush;
}
//=============================================================================
#endif
//=============================================================================
