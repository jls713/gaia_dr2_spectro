#ifndef MASSFUNCTION_H
#define MASSFUNCTION_H
//=============================================================================
// Mass function
//=============================================================================
// Defines the choice of mass function to use in sampling from the EDF
// Importantly there is a distinction between sampling from the initial mass
// function and the present mass function.
// When sampling from the IMF(M), the EDF describes the distribution
// of all stars assuming infinite lifetimes. The corresponding star-formation
// history function Gamma(tau) is the star formation rate.
// When sampling from the PMF(M,tau) (that is the
// IMF/integral(IMF, Min Mass, Max Mass given age) ) the EDF describes the
// distribution of all stars today (excluding those that have died). This means
// the corresponding 'sfh' is not the star-formation history but the age
// distribution of stars that still exist today. To find the SFH we must
// multiply by integral(IMF, Min Mass, Max Mass given age).
//
// To find Max mass at age (also weak function of metallicity), we use
// functionality in chem_evo
//=============================================================================
#include "chemevo/inc/params.h"
#include "chemevo/inc/imf.h"
#include "chemevo/inc/ages.h"
//=============================================================================
class MassFunction{
private:
	bool IMFflag;									// if true, use IMF not PMF
	std::shared_ptr<InitialMassFunction> IMF;
	std::shared_ptr<StellarLifetime> Lifetime;
	const unsigned NM = 1000;
	VecDoub LogMassGrid, CMFgrid;
	double Zsun = 0.017;
public:
	MassFunction(json parameters);
	double CMF(double MaxMass);
	double MF(double M, double tau, double m_h){
		double Z = pow(10.,m_h)*Zsun;
		if(IMFflag)
			return (*IMF)(M);
		else{
			double Mdie = Lifetime->mass_star_dying_now(tau, Z);
			if(M>Mdie)
				return 0.;
			return (*IMF)(M)/ CMF(Mdie);
		}
	}
};
#endif
