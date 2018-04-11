#include "iarates.h"
//=============================================================================
double _binary_sfr(double mu, void *p){
    TypeIaRate_BinaryMass_st *P = (TypeIaRate_BinaryMass_st *) p;
    double m2 = mu*P->m;
    return P->tia->f(mu)*P->tia->SFR(P->R,P->t-MIN(P->t,P->tia->lifetime(m2)));
}

double _typeiarate_binary_sfr(double m, void *p){
    TypeIaRate_BinaryMass_st *P = (TypeIaRate_BinaryMass_st *) p;
    GaussLegendreIntegrator GL(50);
    TypeIaRate_BinaryMass_st PP={P->tia,P->R,P->t,m};
    double mu_b_min=P->tia->min_mass_fraction(P->t,m);
    return P->tia->IMF(m)*GL.integrate(&_binary_sfr,mu_b_min,0.5,&PP);
}
double _binary_RM_sfr(double mu, void *p){
    TypeIaRate_BinaryMass_RM_st *P = (TypeIaRate_BinaryMass_RM_st *) p;
    double m2 = mu*P->m;
    return P->tia->f(mu)*P->tia->SFR(P->R,P->t-MIN(P->t,P->tia->lifetime(m2)))*P->tia->RMKernel(P->R,P->Rp,MIN(P->t,P->tia->lifetime(m2)));
}
int _typeiarate_binary_sfr_2D(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
    double y2[2];
    TypeIaRate_BinaryMass_st_2D *P = (TypeIaRate_BinaryMass_st_2D *) fdata;
    for(int i=0;i<2;i++)
        y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
    auto Mp = y2[0], Rp = y2[1];
    GaussLegendreIntegrator GL(50);
    TypeIaRate_BinaryMass_RM_st PP={P->tia,Rp,P->t,Mp,Rp};
    double mu_b_min=P->tia->min_mass_fraction(P->t,Mp);
    fval[0]=P->tia->IMF(Mp)*GL.integrate(&_binary_RM_sfr,mu_b_min,0.5,&PP);
    return 0;
}

double TypeIaRate_BinaryMass::min_mass_fraction(double t,double m){
    // second mass must be at most as old as the Galaxy.
    // but also the mass of the primary must be less than 8M_solar = 1/2 MB_max
    auto M2 = stellar_lifetime->mass_star_dying_now(t,0.02);
    double mu_min;
    if(m-M2>0.5*MB_max) mu_min=(m-0.5*MB_max)/m;
    else mu_min=M2/m;
    if(mu_min>0.5) return 0.5;
    else if(mu_min<0.) return 0.;
    else return mu_min;
}
double TypeIaRate_BinaryMass::operator()(double R, double t){
    if(t<tau_min) return 0.;
    if(typeid(*radial_migration)==typeid(NoneRadialMigration)){
        GaussLegendreIntegrator GL(100);
        TypeIaRate_BinaryMass_st T={this,R,t,0};
        return ABinary*GL.integrate(&_typeiarate_binary_sfr,MB_min,MB_max,&T);
    }
    else{
        double Rmin=0.,Rmax=30.;
        if(typeid(*radial_migration)==typeid(GaussianRadialMigration)){
            double sR = std::dynamic_pointer_cast<GaussianRadialMigration>(radial_migration)->sigmaR();
            Rmin=R-4.*sR;
            Rmin=(Rmin<0.?0.:Rmin);
            Rmax=R+4.*sR;
        }
        TypeIaRate_BinaryMass_st_2D P(this,R,t,0,{MB_min,Rmin},{MB_max,Rmax});
        double err;
        return ABinary*integrate(&_typeiarate_binary_sfr_2D,&P,1e-4,0,"Divonne",&err,"TypeIaRate_BinaryMass ");
    }
}

double TypeIaRate_BinaryMass::f(double mu){
    return pow(2.*mu,gamma)*2.*(1.+gamma);
}
//=============================================================================
double _typeia_dtd_integrand(double tau, void *P){
    TypeIaRate_DTD_st *T = (TypeIaRate_DTD_st *)P;
    return T->tia->DTD(tau)*T->tia->SFR(T->R,T->t-tau);
}
int _typeia_dtd_integrand_2D(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
    double y2[2];
    TypeIaRate_DTD_st_2D *P = (TypeIaRate_DTD_st_2D *) fdata;
    for(int i=0;i<2;i++)
        y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
    auto tp = exp(y2[0]), Rp = y2[1];
    fval[0]=P->tia->DTD(tp)*P->tia->SFR(Rp,P->t-tp)*P->tia->RMKernel(P->R,Rp,tp)*tp;
    return 0;
}
double norm(double t,void *P){
    TypeIaRate_DTD *T = (TypeIaRate_DTD *)P;
    return T->DTD(t);
}

double TypeIaRate_DTD::operator()(double R, double t){
    if(t<tau_min) return 0.;
    if(typeid(*radial_migration)==typeid(NoneRadialMigration)){
        GaussLegendreIntegrator GL(100);
        TypeIaRate_DTD_st T={this,R,t};
        return ABinary*GL.integrate(&_typeia_dtd_integrand,tau_min,MIN(tau_max,t),&T);
    }
    else{
        double Rmin=0.,Rmax=30.;
        // if(typeid(*radial_migration)==typeid(GaussianRadialMigration)){
        //     double sR = std::dynamic_pointer_cast<GaussianRadialMigration>(radial_migration)->sigmaR();
        //     Rmin=R-4.*sR;
        //     Rmin=(Rmin<0.?0.:Rmin);
        //     Rmax=R+4.*sR;
        // }
        TypeIaRate_DTD_st_2D P(this,R,t,{log(tau_min),Rmin},{log(MIN(tau_max,t)),Rmax});
        double err;
        return ABinary*integrate(&_typeia_dtd_integrand_2D,&P,1e-4,0,"Divonne",&err,"TypeIaRate_DTD ");
    }
}
//=============================================================================
MVP06_TypeIaRate::MVP06_TypeIaRate(ModelParameters M,std::shared_ptr<InitialMassFunction> imf, std::shared_ptr<StarFormationRate> sfr,std::shared_ptr<RadialMigration> RM, std::shared_ptr<StellarLifetime> SL):TypeIaRate_DTD(M,imf,sfr,RM,SL){
    GaussLegendreIntegrator GL(50);
    Norm=1.;
    Norm = GL.integrate(&norm,tau_min,tau_max,this);
}

double MVP06_TypeIaRate::DTD(double t){
    if(t==0.) return 0.;
    double lt = log10(t);
    if(lt<7.93-9.) return pow(10.,1.4-50.*pow(lt+9.-7.7,2.))/Norm;
    else return pow(10.,-0.8-0.9*pow(lt+9.-8.7,2.))/Norm;
}
//=============================================================================
unique_map< TypeIaRate,ModelParameters,
            std::shared_ptr<InitialMassFunction>,
            std::shared_ptr<StarFormationRate>,
            std::shared_ptr<RadialMigration>,
            std::shared_ptr<StellarLifetime>> tia_types ={
    {"Binary",&createInstance<TypeIaRate,TypeIaRate_BinaryMass>},
    {"Matteucci2006",&createInstance<TypeIaRate,MVP06_TypeIaRate>}
};
//=============================================================================
