#ifndef DISTANCES_H
#define DISTANCES_H
//=============================================================================
#include "iso_grid.h"
#include "prior.h"
//=============================================================================
inline double ND_GFunction(VecDoub means, MatDoub icov, double fac){
    return fac*exp(-0.5*(means*(icov*means)));
}

inline double logGFunction(double y, double s){
    s*=s;
    return -.5*log(2.*PI*s)-.5*y*y/s;
}

inline double logND_GFunction(VecDoub means, MatDoub icov, double fac){
    return log(fac)-0.5*(means*(icov*means));
}
/**
 * @brief Class for computing distances from isochrone grid
 */
template<class T>
class DistanceCalculator{
private:
	isochrone_grid<T> *iso_grid;
public:
//=============================================================================
	DistanceCalculator(isochrone_grid<T> *iso_grid):iso_grid(iso_grid){}
    /**
     * @brief compute ln (prob distance modulus | data) from ZAMS
     * @detail no spectroscopic info and assuming age \approx 0
     *
     * @param lbs (Galactic latitude, longitude and distance)
     * @param Z_mass ('true' metallicity, mass (and extinct if required))
     * @param mag (magnitudes)
     * @param err_mag errors in magnitudes
     * @param mag_bands vector of strings giving bands
     * @param bprior turn prior on/off
     * @param EM extinction map pointer
     * @return ln prob(distance modulus|data)
     */
    double photometric_distance_zero_age_pdf(VecDoub lbs, VecDoub Z_mass, VecDoub mag, VecDoub err_mag, std::vector<std::string> mag_bands, galaxy_prior *prior, std::shared_ptr<extinction_map> EM=nullptr);
    /**
     * @brief compute ln (prob distance modulus | data) from ZAMS
     * @detail no spectroscopic info and assuming age \approx 0
     *
     * @param lbs (Galactic latitude, longitude and distance)
     * @param Z_age_mass ('true' metallicity, age and mass (and extinct if required))
     * @param mag (magnitudes)
     * @param err_mag errors in magnitudes
     * @param mag_bands vector of strings giving bands
     * @param bprior turn prior on/off
     * @param EM extinction map pointer
     * @return ln prob(distance modulus|data)
     */
    double photometric_distance_pdf(VecDoub lbs, VecDoub Z_age_mass, VecDoub mag, VecDoub err_mag, std::vector<std::string> mag_bands, galaxy_prior *prior, std::shared_ptr<extinction_map> EM=nullptr);
    /**
     * @brief computes ln prob(distance modulus|data)
     * @details computes the probability of a distance given the spectroscopic
     * data and some estimate of the age, metallicity and mass
     *
     * @param ZTG (measured metallicity, temperature, log g)
     * @param lbs (Galactic latitude, longitude and distance)
     * @param icov inverse covariance matrix of metallicity, teff and logg
     * @param Z_age_mass_alpha ('true' metallicity, age, mass (and alpha))
     * @param mag (magnitudes)
     * @param err_mag errors in magnitudes
     * @param mag_bands vector of strings giving bands
     * @param bprior turn prior on/off
     * @param EM extinction map pointer
     * @return ln prob(distance modulus|data)
     */
    double distance_pdf(VecDoub ZTG, VecDoub lbs, VecDoub icov, VecDoub Z_age_mass_alpha, VecDoub mag, VecDoub err_mag, std::vector<std::string> mag_bands, galaxy_prior *prior, std::shared_ptr<extinction_map> EM=nullptr);
    // as above but with Z_age_mass_alpha=(Z,age,evo stage, alpha)
    double distance_pdf_es(VecDoub ZTG, VecDoub lbs, VecDoub icov, VecDoub Z_age_mass_alpha, VecDoub mag, VecDoub err_mag, std::vector<std::string> mag_bands, galaxy_prior *prior, std::shared_ptr<extinction_map> EM=nullptr);
//=============================================================================
    /**
     * @brief distance estimator using unextincted magnitudes
     * @details computes an estimate of the distance from spectrophotometric
     *          data. Returns a vector of (\int p(s) ds, <DM>, sigma(DM), <s>,
     *          sigma(s), <parallax>, sigma(parallax), <age>, sigma(age),
     *          <mass>, sigma(mass), <metallicity>, sigma(metallicity)).
     *          if return_pdf is true, returns a vector of p(mu) where the
     *          final 8 elements are (mu[0],delta mu, <age>, sigma(age),
     *          <mass>, sigma(mass), <metallicity>, sigma(metallicity)))
     *
     * @param mag vector of magnitudes
     * @param Z metallicity
     * @param Teff log10(Teff)
     * @param logg surface gravity
     * @param l Galactic longitude in radians
     * @param b Galactic latitude in radians
     * @param err_mag vector of errors in magnitudes
     * @param cov covariance matrix in metallicity, Teff and logg
     *          (C_ZZ, C_ZT, C_ZG, C_TT, C_TG, C_GG)
     * @param bprior uses prior if True
     * @param return_pdf returns full pdf (plus extras, see above) if True
     * @param NSTD number of standard deviations away from values to consider
     * @param mag_list vector of strings of magnitude bands
     * @param parallax -- measurement of parallax
     * @param parallax_error -- measurement of parallax error
     *          (if <0 parallax not used)
     * @param mass -- independent measurement of mass (e.g. from [C/N])
     * @param mass_error -- measurement of mass error
     *          (if <0 mass not used)
     * @return vector described above
     */
    VecDoub prob_distance(VecDoub mag,double Z,double Teff,double logg,double l,double b,VecDoub err_mag,VecDoub cov,galaxy_prior *prior,  bool return_pdf, double NSTD,std::vector<std::string> mag_list, double parallax=0., double parallax_error=-1., double mass=0., double mass_error=-1.);
    /**
     * @brief distance estimator using extincted magnitudes
     * @details computes an estimate of the distance from spectrophotometric
     *          data. Returns a vector of (\int p(s) ds, <DM>, sigma(DM), <s>,
     *          sigma(s), <parallax>, sigma(parallax), <age>, sigma(age),
     *          <mass>, sigma(mass), <metallicity>, sigma(metallicity), <Av>,
     *          sigma(Av)).
     *
     * @param mag vector of magnitudes
     * @param Z metallicity
     * @param Teff log10(Teff)
     * @param logg surface gravity
     * @param l Galactic longitude in radians
     * @param b Galactic latitude in radians
     * @param err_mag vector of errors in magnitudes
     * @param covar_ZTG covariance matrix in metallicity, Teff and logg
     *          (C_ZZ, C_ZT, C_ZG, C_TT, C_TG, C_GG)
     * @param bprior uses prior if True
     * @param NSTD number of standard deviations away from values to consider
     * @param mag_list vector of strings of magnitude bands
     * @param EM pointer to extinction map instance
     * @param Aprior initial guess of extinction
     * @param parallax -- measurement of parallax
     * @param parallax_error -- measurement of parallax error
     *          (if <0 parallax not used)
     * @param mass -- independent measurement of mass (e.g. from [C/N])
     * @param mass_error -- measurement of mass error
     *          (if <0 mass not used)
     * @return vector described above
     */
    VecDoub prob_distance_extinct(VecDoub mag,double Z,double Teff,double logg,double l,double b,VecDoub err_mag,VecDoub covar_ZTG,galaxy_prior *prior,  bool return_pdf, double NSTD,std::vector<std::string> mag_list,extinction_map *EM,double Aprior=0., double parallax=0., double parallax_error=-1.,double mass=0.,double mass_error=-1.);
    /**
     * @brief distance estimator using unextincted magnitudes
     * @details computes an estimate of the distance from spectrophotometric
     *          data. Returns a vector of (\int p(s) ds, <DM>, sigma(DM), <s>,
     *          sigma(s), <parallax>, sigma(parallax), <age>, sigma(age),
     *          <mass>, sigma(mass), <metallicity>, sigma(metallicity),
     *          <alpha>, sigma(alpha)).
     *
     * @param mag vector of magnitudes
     * @param Z metallicity
     * @param Teff log10(Teff)
     * @param logg surface gravity
     * @param l Galactic longitude in radians
     * @param b Galactic latitude in radians
     * @param err_mag vector of errors in magnitudes
     * @param cov covariance matrix in metallicity, Teff, logg and alpha
     *          (C_ZZ, C_ZT, C_ZG, C_ZA, C_TT, C_TG, C_TA, C_GG, C_GA, C_AA)
     * @param bprior uses prior if True
     * @param NSTD number of standard deviations away from values to consider
     * @param mag_list vector of strings of magnitude bands
     * @param parallax -- measurement of parallax
     * @param parallax_error -- measurement of parallax error
     *          (if <0 parallax not used)
     * @param mass -- independent measurement of mass (e.g. from [C/N])
     * @param mass_error -- measurement of mass error
     *          (if <0 mass not used)
     * @return vector described above
     */
    VecDoub prob_distance_alpha(VecDoub mag,double Z,double Teff,double logg,double alpha,double l,double b,VecDoub err_mag,VecDoub covar_ZTLA,galaxy_prior *prior, double NSTD,std::vector<std::string> mag_list, double parallax=0., double parallax_error=-1.,double mass=0., double mass_error=-1.);
    /**
     * @brief distance estimator using extincted magnitudes
    /**
     * @brief distance estimator using extincted magnitudes
     * @details computes an estimate of the distance from spectrophotometric
     *          data. Returns a vector of (\int p(s) ds, <DM>, sigma(DM), <s>,
     *          sigma(s), <parallax>, sigma(parallax), <age>, sigma(age),
     *          <mass>, sigma(mass), <metallicity>, sigma(metallicity),
     *          <alpha>, sigma(alpha)).
     *
     * @param mag vector of magnitudes
     * @param Z metallicity
     * @param Teff log10(Teff)
     * @param logg surface gravity
     * @param l Galactic longitude in radians
     * @param b Galactic latitude in radians
     * @param err_mag vector of errors in magnitudes
     * @param cov covariance matrix in metallicity, Teff, logg and alpha
     *          (C_ZZ, C_ZT, C_ZG, C_ZA, C_TT, C_TG, C_TA, C_GG, C_GA, C_AA)
     * @param bprior uses prior if True
     * @param NSTD number of standard deviations away from values to consider
     * @param mag_list vector of strings of magnitude bands
     * @param Aprior_VEC a vector of the mean values for the prior on log AV
     * @param sigAprior_VEC a vector of the std values for the prior on log AV
     * @param log_dist_AV a vector of the distances at which prior values of log AV
     *          are provided
     * @param parallax -- measurement of parallax
     * @param parallax_error -- measurement of parallax error
     *          (if <0 parallax not used)
     * @param mass -- independent measurement of mass (e.g. from [C/N])
     * @param mass_error -- measurement of mass error
     *          (if <0 mass not used)
     * @return vector described above
     */
     VecDoub prob_distance_extinctprior(VecDoub mag,
                                        double Z,double Teff,double logg,
                                        double l,double b,
                                        VecDoub err_mag,
                                        VecDoub covar_ZTG,
                                        galaxy_prior *prior,
                                        bool return_pdf,
                                        double NSTD,
                                        std::vector<std::string> mag_list,
                                        VecDoub Aprior_VEC,
                                        VecDoub sigAprior_VEC,
                                        VecDoub log_dist_AV,
                                        extinction_law *EL,
                                        double parallax=0.,
                                        double parallax_error=-1.,
                                        double mass=0.,
                                        double mass_error=-1.);
};
#endif
