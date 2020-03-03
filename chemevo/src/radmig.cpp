#include "radmig.h"
//=============================================================================
double _convolve(double Rp, void *P){
	convolve_struct *p = (convolve_struct*)P;
	// we use the log extrapolation in grid (second true, first is
	// log interpolation).
	// added the R'/R factor so Sigma(R)2piR = int dR' 2piR' Sigma(R')
	return (Rp/p->R)*(*p->gm)(Rp,p->t,true,true)*(*p->rm)(p->R,Rp,p->dt);
}
double RadialMigration::convolve(Grid*gas_mass, unsigned nR, unsigned nt){
	integrator GL(200000);
	double R = gas_mass->grid_radial()[nR];
	double t = gas_mass->grid_time()[nt-1];
	double dt= gas_mass->grid_time()[nt]-t;
	convolve_struct P(this,gas_mass,t,R,dt);
	double IE=1e-4;

	double Ru = gas_mass->grid_radial().back();
	double Rd = gas_mass->grid_radial().front();
	Rd=0.005;Ru+=3.;

	double interval=1., Rd1=R-interval,Ru1=R+interval;
	if(Ru1>Ru) Ru1=Ru;
	if(Rd1<Rd) Rd1=Rd;
	double integral = GL.integrate(&_convolve,Rd1,Ru1,IE,&P);
	double extra_integral=integral;
	while(fabs(extra_integral)>IE*fabs(integral)){
		interval*=2.;
		double Rd2=R-interval,Ru2=R+interval;
		if(Ru2>Ru) Ru2=Ru;
		if(Rd2<Rd) Rd2=Rd;
		extra_integral=GL.integrate(&_convolve,Rd2,Rd1,IE,&P)
			      +GL.integrate(&_convolve,Ru1,Ru2,IE,&P);
		Rd1=Rd2;Ru1=Ru2;
		integral+=extra_integral;
	}
	// double integral = GL.integrate(&_convolve,Rd,Ru,IE,&P);
	return integral;
}
double _convolve_mass_frac(double Rp, void *P){
	convolve_mf_struct *p = (convolve_mf_struct*)P;
	// we use the log extrapolation in grid (second true, first is
	// log interpolation).
	return (*p->gm)(Rp,p->t,true,true)*(*p->mf)(Rp,p->t,false,false)*(*p->rm)(p->R,Rp,p->dt);
}
double RadialMigration::convolve_massfrac(Grid*gas_mass, Grid*mass_fraction, unsigned nR, unsigned nt){
	integrator GL(200000);
	double R = gas_mass->grid_radial()[nR];
	double t = gas_mass->grid_time()[nt-1];
	double dt= gas_mass->grid_time()[nt]-t;
	convolve_mf_struct P(this,gas_mass,mass_fraction,t,R,dt);
	double IE=1e-4;
	double Ru = gas_mass->grid_radial().back();
	double Rd = gas_mass->grid_radial().front();
	Rd=0.005;Ru+=3.;
	//return GL.integrate(&_convolve_mass_frac,Rd,R,1e-3,&P)+GL.integrate(&_convolve_mass_frac,R,Ru,1e-3,&P);
	double interval=1., Rd1=R-interval,Ru1=R+interval;
	if(Ru1>Ru) Ru1=Ru;
	if(Rd1<Rd) Rd1=Rd;
	double integral = GL.integrate(&_convolve_mass_frac,Rd1,Ru1,IE,&P);
	double extra_integral=integral;
	while(fabs(extra_integral)>IE*fabs(integral)){
		interval*=2.;
		double Rd2=R-interval,Ru2=R+interval;
		if(Ru2>Ru) Ru2=Ru;
		if(Rd2<Rd) Rd2=Rd;
		extra_integral=GL.integrate(&_convolve_mass_frac,Rd2,Rd1,IE,&P)
			      +GL.integrate(&_convolve_mass_frac,Ru1,Ru2,IE,&P);
		Rd1=Rd2;Ru1=Ru2;
		integral+=extra_integral;
	}

	// double integral = GL.integrate(&_convolve_mass_frac,Rd,Ru,IE,&P);
	return integral;
}
//=============================================================================
GaussianRadialMigration::GaussianRadialMigration(ModelParameters M){
	sigmaR0 = extract_param(M.parameters["migration"],"sigmaR",3.);
}

double GaussianRadialMigration::operator()(double R, double Rp, double t){
	auto sigmaR = sigmaR0*sqrt(t);
	sigmaR*=sigmaR;
	double Norm = 2./(1.+erf(Rp/sqrt(2.*sigmaR)));
	return exp(-pow(R-Rp,2.)/2./sigmaR)/sqrt(2.*PI*sigmaR)*Norm;
}
//=============================================================================
GaussianRadialMigration_Drift::GaussianRadialMigration_Drift(ModelParameters M){
    double galaxy_age = M.parameters["fundamentals"]["GalaxyAge"];
	sigmaR0 = sigmaR0 = extract_param(M.parameters["migration"],"sigmaR",3.)/sqrt(galaxy_age);
	Rd = M.parameters["fundamentals"]["StarScaleLength"];
	drift_term = -sigmaR0*sigmaR0*.5/Rd;
}
double GaussianRadialMigration_Drift::operator()(double R, double Rp, double t){
	auto sigmaR = sigmaR0*sqrt(t);
	auto dv = drift_term*t;
	sigmaR*=sigmaR;
	double Norm = 2./(1.+erf((R+dv)/sqrt(2.*sigmaR)));
	return exp(-pow(R-Rp-dv,2.)/2./sigmaR)/sqrt(2.*PI*sigmaR)*Norm;
}
//=============================================================================
// Map for creating shared pointer instances of RadialMigration
// from string of class name
shared_map< RadialMigration,ModelParameters> rm_types ={
    {"Gaussian",&createSharedInstance<RadialMigration,GaussianRadialMigration>},
    {"GaussianDrift",&createSharedInstance<RadialMigration,GaussianRadialMigration_Drift>},
    {"None",&createSharedInstance<RadialMigration,NoneRadialMigration>}
};
//=============================================================================
void CrankNicolsonSolver::step(Grid*gas_mass, unsigned nt, double dt){

	#pragma omp parallel for schedule(dynamic)
	for(unsigned nR=0u;nR<NR;++nR){
		double width, width_down, width_up, gmprev_d, gm_prev,rmlower,rmupper,rmlower_1,rmupper_1,R,Rdown,Rup;

		R = gas_mass->grid_radial()[nR];
		Rdown=gas_mass->Rdown(nR);
		Rup=gas_mass->Rup(nR);

		width_up = gas_mass->annulus_width(nR+1);
		width = gas_mass->annulus_width(nR);
		width_down = gas_mass->annulus_width(nR-1);

		rmlower=(*rad_mig)(Rdown,R,dt)*width*.5; // leaving gas moving in
		rmupper=(*rad_mig)(Rup,R,dt)*width*.5; // leaving gas moving out
		rmlower_1=(*rad_mig)(R,Rdown,dt)*width_down*.5; // entering gas moving in
		rmupper_1=(*rad_mig)(R,Rup,dt)*width_up*.5;

		if(nR==0) rmlower=0.;
		gsl_vector_set(diag_inv,nR,1.+rmupper+rmlower);
		if(nR>0)
			gsl_vector_set(diag_down_inv,nR-1,-rmlower_1);
		if(nR<NR)
			gsl_vector_set(diag_up_inv,nR,-rmupper_1);

		gm_prev = (*gas_mass)(nR,nt-1);
		auto outer_gas_mass=0.;
		if(nR<NR-1)
			outer_gas_mass=(*gas_mass)(nR+1,nt-1);
		else
			outer_gas_mass=gas_mass->log_extrapolate_high(Rup,nt-1);

		if(nR>0){
			gmprev_d=(*gas_mass)(nR-1,nt-1);
			gsl_vector_set(rhs,nR,gmprev_d*rmlower_1+outer_gas_mass*rmupper_1+gm_prev*(1.-rmupper-rmlower));
		}
		else{
			// zero flux boundary condition
			gsl_vector_set(rhs,0,outer_gas_mass*rmupper_1+gm_prev*(1.-rmupper));
		}
	}
	gsl_linalg_solve_tridiag (diag_inv,diag_up_inv,diag_down_inv,rhs,sol);
	for(unsigned nR=0u;nR<NR;++nR)
		gas_mass->set(gsl_vector_get(sol,nR),nR,nt);
	gas_mass->set(gsl_vector_get(rhs,NR-1)*2.-(*gas_mass)(NR-1,nt-1),NR-1,nt);
}

void ForwardSolver::step(Grid*gas_mass, unsigned nt, double dt){
	#pragma omp parallel for schedule(dynamic)
	for(auto nR=0u;nR<NR;++nR){

		double width, width_down, width_up, gmprev_d, gm, gm_prev,rmlower,rmupper,rmlower_1,rmupper_1,R,Rdown,Rup;

		R = gas_mass->grid_radial()[nR];
		Rdown=gas_mass->Rdown(nR);
		Rup=gas_mass->Rup(nR);

		width_up = gas_mass->annulus_width(nR+1);
		width = gas_mass->annulus_width(nR);
		width_down = gas_mass->annulus_width(nR-1);

		gm = (*gas_mass)(nR,nt-1);
		gm_prev = (*gas_mass)(nR,nt-1);

		auto outer_gas_mass=0.;
		if(nR<NR-1)
			outer_gas_mass=(*gas_mass)(nR+1,nt-1);
		else
			outer_gas_mass=gas_mass->log_extrapolate_high(Rup,nt-1);

		width_up = gas_mass->annulus_width(nR+1);
		width = gas_mass->annulus_width(nR);
		width_down = gas_mass->annulus_width(nR-1);
		rmlower=(*rad_mig)(Rdown,R,dt)*width; // leaving gas moving in
		rmupper=(*rad_mig)(Rup,R,dt)*width; // leaving gas moving out
		rmlower_1=(*rad_mig)(R,Rdown,dt)*width_down; // entering gas moving in
		rmupper_1=(*rad_mig)(R,Rup,dt)*width_up; // entering gas moving out
		if(nR>0){
			gmprev_d=(*gas_mass)(nR-1,nt-1);
			gm+=gmprev_d*rmlower_1+outer_gas_mass*rmupper_1-gm_prev*(rmupper+rmlower);
		}
		else{
			// zero flux boundary condition
			gm+=outer_gas_mass*rmupper_1-gm_prev*rmupper;
		}
		gas_mass->set(gm,nR,nt);
	}
}
