#ifndef RADMIG_H
#define RADMIG_H
//=============================================================================
#include "utils.h"
#include "params.h"
#include "grid.h"
#include "gsl/gsl_linalg.h"
//=============================================================================
/**
 * @brief Base class for radial migration
 *
 * dLz f(Lz) = dR R \Sigma(R) so there is a factor of (Rp/R) in the kernel
 */
class RadialMigration{
public:
	/**
	 * @brief Radial migration kernel
	 * @details probability of star migrating from radius Rp to radius R in
	 * time t
	 *
	 * @param R radius now
	 * @param Rp birth radius
	 * @param t time
	 * @return radial migration probability
	 */
	virtual double operator()(double R, double Rp, double t)=0;
	/**
	 * @brief convolve grid at time nt-1 to find value at nR grid point at time nt
	 * @param gm grid (gas mass)
	 * @param nR radial grid point to compute convolution
	 * @param nt time grid point for result
	 * @param dt time interval over which to compute -- defaults to one time step
	 * @return convolution result
	 */
	double convolve(Grid* gm, unsigned nR, unsigned nt, double dt=-1);
	/**
	 * @brief convolve grid at time nt-1 to find value at nR grid point at time nt
	 * @param gm grid (gas mass)
	 * @param mf mass fraction for abundance
	 * @param nR radial grid point to compute convolution
	 * @param nt time grid point for result
	 * @param dt time interval over which to compute -- defaults to one time step
	 * @return convolution result
	 */
	double convolve_massfrac(Grid* gm, Grid*mf, unsigned nR, unsigned nt, double dt=-1);
};
struct convolve_struct{
	RadialMigration *rm;
	Grid *gm;
	double t, R, dt;
	convolve_struct(RadialMigration *rm, Grid *gm, double t, double R, double dt): rm(rm),gm(gm),t(t),R(R),dt(dt){}
};
struct convolve_mf_struct: convolve_struct{
	Grid*mf;
	convolve_mf_struct(RadialMigration *rm, Grid *gm, Grid *mf, double t, double R, double dt): convolve_struct(rm,gm,t,R,dt),mf(mf){}
};
//=============================================================================
/**
 * @brief No Radial migration
 */
class NoneRadialMigration:public RadialMigration{
public:
	NoneRadialMigration(ModelParameters M){}
	double operator()(double R, double Rp, double t){
		throw std::invalid_argument("None Radial Migration can't be used\n");
	}
};
//=============================================================================
/**
 * @brief Simple Gaussian radial migration kernel
 */
class GaussianRadialMigration:public RadialMigration{
private:
	double sigmaR0; // Radial width of stars born at single radius after 1 Gyr
public:
	GaussianRadialMigration(ModelParameters M);
	double sigmaR(void){return sigmaR0;}
	double operator()(double R, double Rp, double t);
};
//=============================================================================
/**
 * @brief Gaussian radial migration kernel with drift to preserve exp profile
 */
class GaussianRadialMigration_Drift:public RadialMigration{
private:
	double sigmaR0; // Radial width of stars born at single radius after 1 Gyr
	double Rd;		// Scale-length of population
	double drift_term; // = -.5*sigma_R0^2/Rd
public:
	GaussianRadialMigration_Drift(ModelParameters M);
	double sigmaR(void){return sigmaR0;}
	double operator()(double R, double Rp, double t);
};
//=============================================================================
// Map for creating new instances of RadialMigration from a string giving the
// class name rm_types[std::string s="<ClassName>"](ModelParameters M)
// produces a shared pointer of Class <ClassName> initialized with the model
// parameters M.
extern shared_map<RadialMigration,ModelParameters> rm_types;
//=============================================================================
/**
 * @brief Crank Nicolson solver for radial migration
 * @details steps forward the diffusion equation controlling the strength of
 * 			radial migration.
 */
class CrankNicolsonSolver{
private:
	std::shared_ptr<RadialMigration> rad_mig;		// rad mig kernel
	unsigned NR; // spatial grid size
	gsl_vector *diag_inv,*diag_up_inv,*diag_down_inv;// stores tridiagonal elements of the propagation matrices
	gsl_vector *rhs, *sol; // rhs = P x (where P is the half-propagation matrix), sol = P'^{-1} P x where P' is the other half-propagation matrix
public:
	/**
	 * @brief Crank Nicolson solver
	 *
	 * @param rad_mig radial migration kernel
	 * @param NR size of spatial grid
	 */
	CrankNicolsonSolver(std::shared_ptr<RadialMigration> rad_mig,unsigned NR):
		rad_mig(rad_mig),NR(NR){
		diag_inv = gsl_vector_alloc(NR);
		diag_up_inv = gsl_vector_alloc(NR-1);
		diag_down_inv = gsl_vector_alloc(NR-1);
		rhs = gsl_vector_alloc(NR);
		sol = gsl_vector_alloc(NR);
	}
	// Destructor
	~CrankNicolsonSolver(){
		gsl_vector_free (diag_inv);
		gsl_vector_free (diag_up_inv);
		gsl_vector_free (diag_down_inv);
		gsl_vector_free (rhs);
		gsl_vector_free (sol);
	}
	/**
	 * @brief step forward C-N solver by one step of dt from time grid point nt to nt+1
	 * @param gas_mass gas mass grid
	 * @param nt current time element
	 * @param dt time step
	 */
	void step(Grid* gas_mass, unsigned nt, double dt);
};
//=============================================================================
/**
 * @brief Forward integration solver for radial migration
 * @details steps forward the diffusion equation controlling the strength of
 * 			radial migration.
 */
class ForwardSolver{
private:
	std::shared_ptr<RadialMigration> rad_mig;		// rad mig kernel
	unsigned NR; // spatial grid size
public:
	/**
	 * @brief Crank Nicolson solver
	 *
	 * @param rad_mig radial migration kernel
	 * @param NR size of spatial grid
	 */
	ForwardSolver(std::shared_ptr<RadialMigration> rad_mig,unsigned NR):
		rad_mig(rad_mig),NR(NR){}
	/**
	 * @brief step forward solver by one step of dt from time grid point nt to nt+1
	 * @param gas_mass gas mass grid
	 * @param nt current time element
	 * @param dt time step
	 */
	void step(Grid* gas_mass, unsigned nt, double dt);
};

#endif
//=============================================================================

