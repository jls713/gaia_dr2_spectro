#include "distances.h"

double salaris_alpha(double a){
    // metallicity difference of isochrone with alpha-enhancement
    return log10(0.638*pow(10.,a)+0.362);
}

//=============================================================================
template<class isochrone_g>
double DistanceCalculator<isochrone_g>::photometric_distance_zero_age_pdf(VecDoub lbs, VecDoub Z_mass, VecDoub mag, VecDoub err_mag, std::vector<std::string> mag_bands, galaxy_prior *prior, std::shared_ptr<extinction_map> EM){

    double Mag, result=0.;
    double mu = 5.*log10(lbs[2]*100.);

    bool extinct = (Z_mass.size()==3);

    double av;
    if(extinct)
        av=exp(Z_mass[2]);
    double metal = Z_mass[0];

    VecDoub iso_r = iso_grid->interp_es(metal,
                                     0.0001,
                                     Z_mass[1],
                                     mag_bands);

    if(iso_r[0]<0.)
        return -std::numeric_limits<double>::infinity();
    if(prior->bprior()){
        VecDoub X = conv::GalacticToCartesian(lbs);
        result+=log(prior->prior(X,
                           metal,
                           0.0001));
    }
    double mag_ex=0.;
    for(unsigned c = 0; c<mag_bands.size(); c++){
        if(extinct) mag_ex=av*EM->extinct_const(mag_bands[c]);
        result+=logGFunction(mag[c]-(mu+iso_r[2+c]+mag_ex),err_mag[c]);
    }

    result+=log(KroupaIMF_default(Z_mass[1]));

    if(EM and extinct)
        result+=log(GFunction(log(av/EM->A_V(lbs[0],lbs[1],lbs[2])),1.)/av);

    double jac = lbs[2];jac*=jac*jac;
    result+=log(jac);
    if(result!=result or std::isnan(result))
        return -std::numeric_limits<double>::infinity();
    return result;
}

//=============================================================================
template<class isochrone_g>
double DistanceCalculator<isochrone_g>::photometric_distance_pdf(VecDoub lbs, VecDoub Z_age_mass, VecDoub mag, VecDoub err_mag, std::vector<std::string> mag_bands, galaxy_prior *prior, std::shared_ptr<extinction_map> EM){

    double Mag, result=0.;
    double mu = 5.*log10(lbs[2]*100.);

    bool extinct = (Z_age_mass.size()==4);

    double av;
    if(extinct)
        av=exp(Z_age_mass[3]);
    double metal = Z_age_mass[0];

    VecDoub iso_r = iso_grid->interp_es(metal,
                                     Z_age_mass[1],
                                     Z_age_mass[2],
                                     mag_bands);

    if(iso_r[0]<0.)
        return -std::numeric_limits<double>::infinity();
    if(prior->bprior()){
        VecDoub X = conv::GalacticToCartesian(lbs);
        result+=log(prior->prior(X,
                           metal,
                           Z_age_mass[1]));
    }
    double mag_ex=0.;
    for(unsigned c = 0; c<mag_bands.size(); c++){
        if(extinct) mag_ex=av*EM->extinct_const(mag_bands[c]);
        result+=logGFunction(mag[c]-(mu+iso_r[2+c]+mag_ex),err_mag[c]);
    }

    if(EM and extinct)
        result+=log(GFunction(log(av/EM->A_V(lbs[0],lbs[1],lbs[2])),1.)/av);

    result+=log(KroupaIMF_default(Z_age_mass[2]));

    double jac = lbs[2];jac*=jac*jac;
    result+=log(jac);
    if(result!=result or std::isnan(result))
        return -std::numeric_limits<double>::infinity();
    return result;
}

//=============================================================================

template<class isochrone_g>
double DistanceCalculator<isochrone_g>::distance_pdf(VecDoub ZTG, VecDoub lbs, VecDoub icov, VecDoub Z_age_mass_alpha, VecDoub mag, VecDoub err_mag, std::vector<std::string> mag_bands, galaxy_prior *prior, std::shared_ptr<extinction_map> EM){

    double Mag, result=0.;
    double mu = 5.*log10(lbs[2]*100.);
    bool alpha=(ZTG.size()==4);
    if(alpha and Z_age_mass_alpha.size()<4)
	throw std::invalid_argument("alpha val passed but no true alpha");
    MatDoub icov_mat(3+alpha,VecDoub(3+alpha,0.));
    icov_mat[0][0]=icov[0];
    icov_mat[0][1]=icov[1];
    icov_mat[1][0]=icov[1];
    icov_mat[0][2]=icov[2];
    icov_mat[2][0]=icov[2];
    icov_mat[1][1]=icov[3+alpha];
    icov_mat[1][2]=icov[4+alpha];
    icov_mat[2][1]=icov[4+alpha];
    icov_mat[2][2]=icov[5+alpha*2];
    if(alpha){
    	icov_mat[3][0]=icov[3];
    	icov_mat[0][3]=icov_mat[3][0];
        icov_mat[3][1]=icov[6];
        icov_mat[1][3]=icov_mat[3][1];
        icov_mat[3][2]=icov[8];
        icov_mat[2][3]=icov_mat[3][2];
    	icov_mat[3][3]=icov[9];
    }

    bool extinct = (Z_age_mass_alpha.size()==(4+alpha));

    double av;
    if(extinct)
        av=exp(Z_age_mass_alpha[3]);
    double metal = Z_age_mass_alpha[0];
    if(alpha) metal+= salaris_alpha(Z_age_mass_alpha[3+extinct]);

    VecDoub iso_r = iso_grid->interp_es(metal,
                                     Z_age_mass_alpha[1],
                                     Z_age_mass_alpha[2],
                                     mag_bands);


    if(iso_r[0]<0.)
        return -std::numeric_limits<double>::infinity();
    if(prior->bprior()){
        VecDoub X = conv::GalacticToCartesian(lbs);
        if(alpha)
    		result+=log(prior->prior(X,
                               Z_age_mass_alpha[0],
                               Z_age_mass_alpha[1],
    			   Z_age_mass_alpha[3+extinct]));
        else
        	result+=log(prior->prior(X,
                           Z_age_mass_alpha[0],
                           Z_age_mass_alpha[1]));
    }
    double mag_ex=0.;
    for(unsigned c = 0; c<mag_bands.size(); c++){
        if(extinct) mag_ex=av*EM->extinct_const(mag_bands[c]);
        result+=logGFunction(mag[c]-(mu+iso_r[2+c]+mag_ex),err_mag[c]);
    }
    VecDoub mshift=  {Z_age_mass_alpha[0]-ZTG[0],
                      iso_r[0]-ZTG[1],
                      iso_r[1]-ZTG[2]};
    if(alpha) mshift.push_back(Z_age_mass_alpha[3+extinct]-ZTG[3]);

    result+=log(KroupaIMF_default(Z_age_mass_alpha[2]))+logND_GFunction(mshift,icov_mat,1.);

    if(EM and extinct)
        result+=log(GFunction(log(av/EM->A_V(lbs[0],lbs[1],lbs[2])),1.)/av);

    double jac = lbs[2];jac*=jac*jac;
    result+=log(jac);
    if(result!=result or std::isnan(result))
        return -std::numeric_limits<double>::infinity();
    return result;
}

//=============================================================================
template<class isochrone_g>
double DistanceCalculator<isochrone_g>::distance_pdf_es(VecDoub ZTG, VecDoub lbs, VecDoub icov, VecDoub Z_age_mass_alpha, VecDoub mag, VecDoub err_mag, std::vector<std::string> mag_bands, galaxy_prior *prior, std::shared_ptr<extinction_map> EM){
    double desdM = iso_grid->d_es_d_M(Z_age_mass_alpha[0],Z_age_mass_alpha[1],Z_age_mass_alpha[2]);
    Z_age_mass_alpha[2] = iso_grid->mass_from_es(Z_age_mass_alpha[0],Z_age_mass_alpha[1],Z_age_mass_alpha[2]);
    return distance_pdf(ZTG, lbs, icov, Z_age_mass_alpha, mag, err_mag, mag_bands, prior, EM)-log(desdM);
}

VecDoub compute_mean_std_from_total(VecDoub total){
    for(auto i=1;i<total.size();++i) total[i]/=total[0];
    for(auto i=2;i<total.size();i+=2)
        total[i]=sqrt(total[i]-total[i-1]*total[i-1]);
    return total;
}

//=============================================================================

template<class isochrone_g>
VecDoub DistanceCalculator<isochrone_g>::prob_distance(
    VecDoub mag,
    double Z,
    double Teff,
    double logg,
    double l,
    double b,
    VecDoub err_mag,
    VecDoub covar_ZTL,
    galaxy_prior *prior,
    bool return_pdf,
    double NSTD,
    std::vector<std::string> mag_list,
    double parallax, double parallax_error,
    double mass, double mass_error){

    double dmu=err_mag[0]/5.,s,Teffmodel,loggmodel,agemodel,mmodel;
    MatDoub icov(3,VecDoub(3,0.)),cov(3,VecDoub(3,0.));
    cov[0][0]=covar_ZTL[0];
    cov[0][1]=covar_ZTL[1]; cov[1][0]=covar_ZTL[1];
    cov[0][2]=covar_ZTL[2]; cov[2][0]=covar_ZTL[2];
    cov[1][1]=covar_ZTL[3];
    cov[1][2]=covar_ZTL[4]; cov[2][1]=covar_ZTL[4];
    cov[2][2]=covar_ZTL[5];
    icov = inverse3D(cov);
    double det = CalcDeterminant(cov,3);
    double fac = det;
    for(unsigned i=0;i<3;++i) fac*=TPI;
    fac = 1./sqrt(fac);
    double err_Z = sqrt(cov[0][0]);
    double err_Teff = sqrt(cov[1][1]);
    double err_logg = sqrt(cov[2][2]);

    double Mag;

    double nstd=NSTD;
    VecDoub total(17,0.); double tmp1, tmp2, tmp3, tmp4, s2;
    VecDoub prob_mu,mu_grid;
    double grid_space=0.2;

    for(int i=0;i<iso_grid->NF;i++){
        if(fabs(iso_grid->iso(i,0)->feh()-Z)>nstd*err_Z) continue;
        tmp1 = iso_grid->delta_Z(i);
        for(int j=0;j<iso_grid->NA;j++){
            agemodel=iso_grid->iso(i,j)->tau();
            tmp2 = tmp1*iso_grid->delta_age(j);
            for(int k=0;k<iso_grid->iso(i,j)->N();k++){

                Teffmodel = iso_grid->iso(i,j)->logTeff(k);
                if(fabs(Teffmodel-Teff)>nstd*err_Teff) continue;
                loggmodel = iso_grid->iso(i,j)->logg(k);
                if(fabs(loggmodel-logg)>nstd*err_logg) continue;

                mmodel=iso_grid->iso(i,j)->initial_mass(k);
                tmp3=tmp2*iso_grid->iso(i,j)->delta_mass(k)
                    *KroupaIMF_default(mmodel)*ND_GFunction({Z-iso_grid->fehgrid[i],Teff-Teffmodel,logg-loggmodel},icov,fac);

                if(mass_error>0.)
                    tmp3 *= GFunction(mass-iso_grid->iso(i,j)->mass(k),
                                      mass_error);

                double MagN=iso_grid->iso(i,j)->mag(k,mag_list[0]);
                for(double  mu=mag[0]-MagN-nstd*err_mag[0];
                            mu<mag[0]-MagN+nstd*err_mag[0];
                            mu+=dmu){
                    s = pow(10.,0.2*mu-2.);s2=s*s;
                    tmp4=s2*s*tmp3*dmu; // p(mu) dmu
                    if(parallax_error>0.)
                        tmp4*=GFunction(parallax-1./s,parallax_error);
                    for(unsigned c = 0; c<mag_list.size(); c++){
                        Mag=iso_grid->iso(i,j)->mag(k,mag_list[c]);
                        tmp4*=GFunction(mag[c]-(mu+Mag),err_mag[c]);
                    }
                    if(tmp4==0.) continue;
                    if(prior->bprior()){
                        VecDoub X = conv::GalacticToCartesian({l,b,s});
                        tmp4*=prior->prior(X,
                                           iso_grid->fehgrid[i],
                                           agemodel);
                    }
                    // Distance modulus
                    total[0]+=tmp4;total[1]+=mu*tmp4;total[2]+=mu*mu*tmp4;
                    // Distance
                    total[3]+=s*tmp4;total[4]+=s2*tmp4;
                    // Parallax
                    total[5]+=tmp4/s;total[6]+=tmp4/s2;
                    // Age
                    total[7]+=tmp4*agemodel;total[8]+=tmp4*agemodel*agemodel;
                    // Mass
                    total[9]+=tmp4*mmodel;total[10]+=tmp4*mmodel*mmodel;
                    // Metallicity
                    total[11]+=tmp4*iso_grid->fehgrid[i];
                    total[12]+=tmp4*iso_grid->fehgrid[i]*iso_grid->fehgrid[i];
                    // Teff
                    total[13]+=tmp4*Teffmodel;total[14]+=tmp4*Teffmodel*Teffmodel;
                    // logg
                    total[15]+=tmp4*loggmodel;total[16]+=tmp4*loggmodel*loggmodel;
                    if(return_pdf){
                        if(prob_mu.size()==0){
                            prob_mu.push_back(tmp4);
                            mu_grid.push_back(mu-grid_space*.5);
                            mu_grid.push_back(mu+grid_space*.5);
                        }
                        else{
                            if(mu<mu_grid.front()){
                                while(mu<mu_grid.front()){
                                    mu_grid.insert(mu_grid.begin(),
                                               mu_grid.front()-grid_space);
                                    prob_mu.insert(prob_mu.begin(),
                                               tmp4);
                                }
                            }
                            else if(mu>mu_grid.back()){
                                while(mu>mu_grid.back()){
                                    mu_grid.push_back(mu_grid.back()+grid_space);
                                    prob_mu.push_back(tmp4);
                                }
                            }
                            else{
                                int bot,top;
                                topbottom(mu_grid,mu,&bot,&top);
                                prob_mu[bot]+=tmp4;
                            }
                        }
                    }
                }
            }
        }
    }
    if(total[0]<=0.){
        if(NSTD>=40.){
            total[0]=std::numeric_limits<double>::quiet_NaN();
            total[1]=std::numeric_limits<double>::quiet_NaN();
            total[2]=-1.;
            for(int i=3;i<13;++i)
                total[i]=std::numeric_limits<double>::infinity();
            return total;
        }
        else return prob_distance(mag,Z,Teff,logg,l,b,err_mag,
                                  covar_ZTL,prior,return_pdf,2*NSTD,
                                  mag_list,parallax,parallax_error);
    }
    total=compute_mean_std_from_total(total);
    if(return_pdf){
        double sum=0.;
        for(auto i:prob_mu)
            sum+=i*grid_space;
        for(unsigned i=0;i<prob_mu.size();++i)
            prob_mu[i]/=sum;
        prob_mu.push_back(mu_grid[0]);
        prob_mu.push_back(mu_grid[1]-mu_grid[0]);
        prob_mu.push_back(total[7]);
        prob_mu.push_back(total[8]);
        prob_mu.push_back(total[9]);
        prob_mu.push_back(total[10]);
        prob_mu.push_back(total[11]);
        prob_mu.push_back(total[12]);
        return prob_mu;
    }
    else
        return total;
}

template<class isochrone_g>
VecDoub DistanceCalculator<isochrone_g>::prob_distance_alpha(
    VecDoub mag,
    double Z,
    double Teff,
    double logg,
    double alpha,
    double l,
    double b,
    VecDoub err_mag,VecDoub covar_ZTLA,
    galaxy_prior *prior,
    double NSTD,
    std::vector<std::string> mag_list,
    double parallax, double parallax_error,
    double mass, double mass_error){

    double dmu=err_mag[0]/5.,s,Teffmodel,loggmodel,agemodel,mmodel,fehmodel;
    MatDoub icov(4,VecDoub(4,0.)),cov(4,VecDoub(4,0.));
    cov[0][0]=covar_ZTLA[0];
    cov[0][1]=covar_ZTLA[1]; cov[1][0]=covar_ZTLA[1];
    cov[0][2]=covar_ZTLA[2]; cov[2][0]=covar_ZTLA[2];
    cov[0][2]=covar_ZTLA[3]; cov[2][0]=covar_ZTLA[3];
    cov[1][1]=covar_ZTLA[4];
    cov[1][2]=covar_ZTLA[5]; cov[2][1]=covar_ZTLA[5];
    cov[1][3]=covar_ZTLA[6]; cov[3][1]=covar_ZTLA[6];
    cov[2][2]=covar_ZTLA[7];
    cov[2][3]=covar_ZTLA[8]; cov[3][2]=covar_ZTLA[8];
    cov[3][3]=covar_ZTLA[9];
    MatrixInversion(cov,icov);
    double det = CalcDeterminant(cov,4);
    double fac = det;
    for(unsigned i=0;i<4;++i) fac*=TPI;
    fac = 1./sqrt(fac);
    double err_Z = sqrt(cov[0][0]);
    double err_Teff = sqrt(cov[1][1]);
    double err_logg = sqrt(cov[2][2]);
    double err_alpha = sqrt(cov[3][3]);

    double Mag;

    double nstd=NSTD;
    VecDoub total(15,0.); double tmp1, tmp2, tmp3, tmp4, s2;
    for(int i=0;i<iso_grid->NF;i++){
        if(fabs(iso_grid->iso(i,0)->feh()-(Z+salaris_alpha(alpha)))>sqrt(2.)*nstd*err_Z) continue;
        tmp1 = iso_grid->delta_Z(i);
        for(int j=0;j<iso_grid->NA;j++){
            agemodel=iso_grid->iso(i,j)->tau();
            tmp2 = tmp1*iso_grid->delta_age(j);
            for(int k=0;k<iso_grid->iso(i,j)->N();k++){

                Teffmodel = iso_grid->iso(i,j)->logTeff(k);
                if(fabs(Teffmodel-Teff)>nstd*err_Teff) continue;
                loggmodel = iso_grid->iso(i,j)->logg(k);
                if(fabs(loggmodel-logg)>nstd*err_logg) continue;

                mmodel=iso_grid->iso(i,j)->initial_mass(k);
                for(double  alphamodel=alpha-5.*err_alpha;
                            alphamodel<alpha+5.*err_alpha;
                            alphamodel+=tmp1){
                    fehmodel = iso_grid->fehgrid[i]-salaris_alpha(alphamodel);
                    tmp3=tmp2*iso_grid->iso(i,j)->delta_mass(k)
                        *KroupaIMF_default(mmodel)
                        *ND_GFunction({Z-fehmodel,
                                       Teff-Teffmodel,
                                       logg-loggmodel,
                                       alpha-alphamodel},icov,fac);
                    if(mass_error>0.)
                        tmp3 *= GFunction(mass-iso_grid->iso(i,j)->mass(k),
                                          mass_error);
                    double MagN=iso_grid->iso(i,j)->mag(k,mag_list[0]);
                    for(double  mu=mag[0]-MagN-nstd*err_mag[0];
                                mu<mag[0]-MagN+nstd*err_mag[0];
                                mu+=dmu){
                        s = pow(10.,0.2*mu-2.);s2=s*s;
                        tmp4=s2*s*tmp3*dmu;  // p(mu) dmu
                        if(parallax_error>0.)
                            tmp4*=GFunction(parallax-1./s,parallax_error);
                        for(unsigned c = 0; c<mag_list.size(); c++){
                            Mag=iso_grid->iso(i,j)->mag(k,mag_list[c]);
                            tmp4*=GFunction(mag[c]-(mu+Mag),err_mag[c]);
                        }
                        if(tmp4==0.) continue;
                        if(prior->bprior()){
                            VecDoub X = conv::GalacticToCartesian({l,b,s});
                            tmp4*=prior->prior(X,
                                               iso_grid->fehgrid[i],
                                               agemodel);
                        }
                        total[0]+=tmp4;total[1]+=mu*tmp4;
                        total[2]+=mu*mu*tmp4;
                        total[3]+=s*tmp4;total[4]+=s2*tmp4;
                        total[5]+=tmp4/s;total[6]+=tmp4/s2;
                        total[7]+=tmp4*agemodel;
                        total[8]+=tmp4*agemodel*agemodel;
                        total[9]+=tmp4*mmodel;
                        total[10]+=tmp4*mmodel*mmodel;
                        total[11]+=tmp4*fehmodel;
                        total[12]+=tmp4*fehmodel*fehmodel;
                        total[13]+=tmp4*alphamodel;
                        total[14]+=tmp4*alphamodel*alphamodel;
                    }
                }
            }
        }
    }
    if(total[0]<=0.){
        if(NSTD>=40.){
            total[0]=std::numeric_limits<double>::quiet_NaN();
            total[1]=std::numeric_limits<double>::quiet_NaN();
            total[2]=-1.;
            for(int i=3;i<15;++i)total[i]=std::numeric_limits<double>::quiet_NaN();return total;
        }
        else return prob_distance_alpha(mag,Z,Teff,logg,alpha,l,b,err_mag,covar_ZTLA,prior,2*NSTD,mag_list,parallax,parallax_error);
    }
    total=compute_mean_std_from_total(total);
    return total;
}


template<class isochrone_g>
VecDoub DistanceCalculator<isochrone_g>::prob_distance_extinct(
    VecDoub mag,
    double Z,
    double Teff,
    double logg,
    double l,
    double b,
    VecDoub err_mag,VecDoub covar_ZTL,
    galaxy_prior *prior,
    bool return_pdf,
    double NSTD,
    std::vector<std::string> mag_list,
    extinction_map *EM,
    double Aprior,
    double parallax, double parallax_error,
    double mass, double mass_error){

    if(Aprior==0.){
        // Only use IR bands
        VecDoub mags_use, err_mags_use; std::vector<std::string> mag_str_use;
        for(unsigned i=0;i<mag_list.size();++i){
            if(mag_list[i]=="H" or mag_list[i]=="J" or mag_list[i]=="K" or mag_list[i]=="r" or mag_list[i]=="i"){
                mags_use.push_back(mag[i]);
                err_mags_use.push_back(err_mag[i]);
                mag_str_use.push_back(mag_list[i]);
            }
        }
        VecDoub first_guess=prob_distance(mags_use,Z,Teff,logg,
                                          l,b,err_mags_use,covar_ZTL,
                                          prior,return_pdf,NSTD,mag_str_use);
        Aprior = EM->A_V(l,b,first_guess[3]);
    }

    VecDoub total(19,0.);
    double dmu=err_mag[0]/2.,s,Teffmodel,loggmodel,agemodel,mmodel;
    MatDoub icov(3,VecDoub(3,0.)),cov(3,VecDoub(3,0.));
    cov[0][0]=covar_ZTL[0];
    cov[0][1]=covar_ZTL[1]; cov[1][0]=covar_ZTL[1];
    cov[0][2]=covar_ZTL[2]; cov[2][0]=covar_ZTL[2];
    cov[1][1]=covar_ZTL[3];
    cov[1][2]=covar_ZTL[4]; cov[2][1]=covar_ZTL[4];
    cov[2][2]=covar_ZTL[5];
    icov = inverse3D(cov);
    double det = CalcDeterminant(cov,3);
    double fac = det;
    for(unsigned i=0;i<3;++i) fac*=TPI;
    fac = 1./sqrt(fac);
    double err_Z = sqrt(cov[0][0]);
    double err_Teff = sqrt(cov[1][1]);
    double err_logg = sqrt(cov[2][2]);
    double Mag_extN;
    double nstd=NSTD;

    double EF;
    VecDoub mag_ext_const(mag_list.size(),0.);
    for(unsigned i=0;i<mag.size();++i)
        mag_ext_const[i]=EM->extinct_const(mag_list[i]);
    double tmp1, tmp2, tmp3, tmp4, s2, Mag, Mag_ext;

    double av = 0.;
    double lAprior = log(Aprior);
    double sigA = (Aprior>0.?0.4/Aprior:0.2);
    double AVmin=(Aprior>0.?lAprior-3.*sigA:0.);
    double AVmax=(Aprior>0.?lAprior+3.*sigA:log(EM->A_V(l,b,100.)));
    if(AVmax>3.4) AVmax=3.4;


    VecDoub prob_mu,mu_grid;
    double grid_space=0.2;

    for(int i=0;i<iso_grid->NF;i++){
        if(fabs(iso_grid->iso(i,0)->feh()-Z)>nstd*err_Z) continue;
        tmp1 = iso_grid->delta_Z(i);
        for(int j=0;j<iso_grid->NA;j++){
            agemodel=iso_grid->iso(i,j)->tau();
            tmp2 = tmp1*iso_grid->delta_age(j);
            for(int k=0;k<iso_grid->iso(i,j)->N();k++){
                Teffmodel = iso_grid->iso(i,j)->logTeff(k);
                if(fabs(Teffmodel-Teff)>nstd*err_Teff) continue;
                loggmodel = iso_grid->iso(i,j)->logg(k);
                if(fabs(loggmodel-logg)>nstd*err_logg) continue;

                mmodel=iso_grid->iso(i,j)->initial_mass(k);
                tmp3=tmp2*iso_grid->iso(i,j)->delta_mass(k)
                    *KroupaIMF_default(mmodel)*ND_GFunction({Z-iso_grid->fehgrid[i],
                                                     Teff-Teffmodel,
                                                     logg-loggmodel},icov,fac);

                if(mass_error>0.)
                    tmp3 *= GFunction(mass-iso_grid->iso(i,j)->mass(k),
                                      mass_error);

                Mag=iso_grid->iso(i,j)->mag(k,mag_list[0]);
                for(double lav=AVmin;lav<AVmax;lav+=sigA/2.){
        		    EF = GFunction(lav-lAprior,sigA);
                    av  = exp(lav);
                    Mag_extN = Mag+av*mag_ext_const[0];
                    for(double mu=mag[0]-Mag_extN-nstd*err_mag[0];
                               mu<mag[0]-Mag_extN+nstd*err_mag[0];
                               mu+=dmu){
                        s = pow(10.,0.2*mu-2.);s2=s*s;
                        tmp4=s2*s*tmp3*dmu*EF;
                        if(parallax_error>0.)
                            tmp4*=GFunction(parallax-1./s,parallax_error);
                        for(unsigned c=0;c<mag.size();++c){
                            Mag_ext=iso_grid->iso(i,j)->mag(k,mag_list[c])+av*mag_ext_const[c];
                            tmp4*=GFunction(mag[c]-(mu+Mag_ext),err_mag[c]);
                        }
                        if(tmp4==0.) continue;
                        if(prior->bprior()){
                            VecDoub X = conv::GalacticToCartesian({l,b,s});
                            tmp4*=prior->prior(X,
                                               iso_grid->fehgrid[i],
                                               agemodel);
                        }
                        total[0]+=tmp4;
                        total[1]+=mu*tmp4;total[2]+=mu*mu*tmp4;
                        total[3]+=s*tmp4;total[4]+=s2*tmp4;
                        total[5]+=tmp4/s;total[6]+=tmp4/s2;
                        total[7]+=tmp4*agemodel;
                        total[8]+=tmp4*agemodel*agemodel;
                        total[9]+=tmp4*mmodel;total[10]+=tmp4*mmodel*mmodel;
                        total[11]+=tmp4*iso_grid->fehgrid[i];total[12]+=tmp4*iso_grid->fehgrid[i]*iso_grid->fehgrid[i];
                        total[13]+=tmp4*av;total[14]+=tmp4*av*av;
                        // Teff
                        total[15]+=tmp4*Teffmodel;
                        total[16]+=tmp4*Teffmodel*Teffmodel;
                        // logg
                        total[17]+=tmp4*loggmodel;
                        total[18]+=tmp4*loggmodel*loggmodel;
                        if(return_pdf){
                            if(prob_mu.size()==0){
                                prob_mu.push_back(tmp4);
                                mu_grid.push_back(mu-grid_space*.5);
                                mu_grid.push_back(mu+grid_space*.5);
                            }
                            else{
                                if(mu<mu_grid.front()){
                                    while(mu<mu_grid.front()){
                                        mu_grid.insert(mu_grid.begin(),
                                                   mu_grid.front()-grid_space);
                                        prob_mu.insert(prob_mu.begin(),
                                                   tmp4);
                                    }
                                }
                                else if(mu>mu_grid.back()){
                                    while(mu>mu_grid.back()){
                                        mu_grid.push_back(mu_grid.back()+grid_space);
                                        prob_mu.push_back(tmp4);
                                    }
                                }
                                else{
                                    int bot,top;
                                    topbottom(mu_grid,mu,&bot,&top);
                                    prob_mu[bot]+=tmp4;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if(total[0]<=0.){
        if(NSTD>=20.){
            total[0]=std::numeric_limits<double>::infinity();
            total[1]=std::numeric_limits<double>::infinity();
        total[2]=-1.;
            for(int i=3;i<13;++i)total[i]=std::numeric_limits<double>::infinity();return total;
    }
        else return prob_distance_extinct(mag,Z,Teff,logg,l,b,err_mag,
                                          covar_ZTL,prior,return_pdf,2*NSTD,
                                          mag_list,EM,Aprior,parallax,
                                          parallax_error);
    }

    total=compute_mean_std_from_total(total);

    if(return_pdf){
        double sum=0.;
        for(auto i:prob_mu)
            sum+=i*grid_space;
        for(unsigned i=0;i<prob_mu.size();++i)
            prob_mu[i]/=sum;
        prob_mu.push_back(mu_grid[0]);
        prob_mu.push_back(mu_grid[1]-mu_grid[0]);
        prob_mu.push_back(total[7]);
        prob_mu.push_back(total[8]);
        prob_mu.push_back(total[9]);
        prob_mu.push_back(total[10]);
        prob_mu.push_back(total[11]);
        prob_mu.push_back(total[12]);
        prob_mu.push_back(total[13]);
        prob_mu.push_back(total[14]);
	return prob_mu;
    }
    else
        return total;
}

template<class isochrone_g>
VecDoub DistanceCalculator<isochrone_g>::prob_distance_extinctprior(
    VecDoub mag,
    double Z,
    double Teff,
    double logg,
    double l,
    double b,
    VecDoub err_mag,VecDoub covar_ZTL,
    galaxy_prior *prior,
    bool return_pdf,
    double NSTD,
    std::vector<std::string> mag_list,
    VecDoub Aprior_VEC,
    VecDoub sigAprior_VEC,
    VecDoub log_dist_AV,
    extinction_law *EL,
    double parallax, double parallax_error,
    double mass, double mass_error){

    double Aprior = 0., sigAprior = 0., lp;

    // Only use IR bands to estimate distance
    VecDoub mags_use, err_mags_use; std::vector<std::string> mag_str_use;
    for(unsigned i=0;i<mag_list.size();++i){
        if(mag_list[i]=="H" or mag_list[i]=="J" or mag_list[i]=="K"
           or mag_list[i]=="r" or mag_list[i]=="i"
           or mag_list[i]=="Jv" or mag_list[i]=="Hv" or mag_list[i]=="Kv"){
            mags_use.push_back(mag[i]);
            err_mags_use.push_back(err_mag[i]);
            mag_str_use.push_back(mag_list[i]);
        }
    }
    if(mag_str_use.size()==0){
        std::cerr<<"No valid bands to initially guess distance.";
        std::cerr<<" Using initial guess of s=1kpc."<<std::endl;
        lp=0.;
    }
    else{
        VecDoub first_guess=prob_distance(mags_use,Z,Teff,logg,l,b,err_mags_use,covar_ZTL,prior,return_pdf,NSTD,mag_str_use,parallax,parallax_error);

        lp = log(first_guess[3]);
    }
    Aprior= linterp(log_dist_AV,Aprior_VEC,lp,"linear");
    sigAprior = linterp(log_dist_AV,sigAprior_VEC,lp,"linear");

    VecDoub total(20,0.);
    double dmu=err_mag[0]/2.,s,Teffmodel,loggmodel,agemodel,mmodel;
    MatDoub icov(3,VecDoub(3,0.)),cov(3,VecDoub(3,0.));
    cov[0][0]=covar_ZTL[0];
    cov[0][1]=covar_ZTL[1]; cov[1][0]=covar_ZTL[1];
    cov[0][2]=covar_ZTL[2]; cov[2][0]=covar_ZTL[2];
    cov[1][1]=covar_ZTL[3];
    cov[1][2]=covar_ZTL[4]; cov[2][1]=covar_ZTL[4];
    cov[2][2]=covar_ZTL[5];
    icov = inverse3D(cov);
    double det = CalcDeterminant(cov,3);
    double fac = det;
    for(unsigned i=0;i<3;++i) fac*=TPI;
    fac = 1./sqrt(fac);
    double err_Z = sqrt(cov[0][0]);
    double err_Teff = sqrt(cov[1][1]);
    double err_logg = sqrt(cov[2][2]);
    double Mag_extN;
    double nstd=NSTD;

    double EF;
    VecDoub mag_ext_const(mag_list.size(),0.);
    // Non-linearity of extinction consts -- only important for G
    VecDoub avgradient(mag_list.size(),0.);
    for(unsigned i=0;i<mag.size();++i){
        mag_ext_const[i]=EL->extinct_const(mag_list[i]);
        avgradient[i]=EL->av_gradient(mag_list[i]);
    }

    double tmp1, tmp2, tmp3, tmp4, s2, Mag, Mag_ext, med_ext, std_ext, DM;

    double av = 0.;
    // double lAprior = log(Aprior);
    // double sigA = (Aprior>0.?0.4/Aprior:0.2);
    double AVmin=Aprior-3.*sigAprior;
    double AVmax=Aprior+3.*sigAprior;
    if(AVmax>Aprior_VEC.back()) AVmax=Aprior_VEC.back();
    if(AVmax>3.4) AVmax=3.4;

    VecDoub prob_mu,mu_grid;
    double grid_space=0.2;
    for(int i=0;i<iso_grid->NF;i++){
        if(fabs(iso_grid->iso(i,0)->feh()-Z)>nstd*err_Z) continue;
        tmp1 = iso_grid->delta_Z(i);
        for(int j=0;j<iso_grid->NA;j++){
            agemodel=iso_grid->iso(i,j)->tau();
            tmp2 = tmp1*iso_grid->delta_age(j);
            for(int k=0;k<iso_grid->iso(i,j)->N();k++){
                Teffmodel = iso_grid->iso(i,j)->logTeff(k);

                for(unsigned i=0;i<mag.size();++i)
                    mag_ext_const[i]=EL->extinct_const(mag_list[i],
                                                       Teffmodel);

                if(fabs(Teffmodel-Teff)>nstd*err_Teff) continue;
                loggmodel = iso_grid->iso(i,j)->logg(k);
                if(fabs(loggmodel-logg)>nstd*err_logg) continue;

                mmodel=iso_grid->iso(i,j)->initial_mass(k);
                tmp3=tmp2*iso_grid->iso(i,j)->delta_mass(k)
                    *KroupaIMF_default(mmodel)*ND_GFunction({Z-iso_grid->fehgrid[i],Teff-Teffmodel,logg-loggmodel},icov,fac);

                if(mass_error>0.)
                    tmp3 *= GFunction(mass-iso_grid->iso(i,j)->mass(k),
                                      mass_error);

                Mag=iso_grid->iso(i,j)->mag(k,mag_list[0]);
                for(double lav=AVmin;lav<AVmax;lav+=sigAprior/2.){
                    av  = exp(lav);
                    Mag_extN = Mag+av*mag_ext_const[0]*(1-av*avgradient[0]);
                    DM = mag[0]-Mag_extN;
                    s = pow(10.,0.2*DM-2.);s2=s*s;lp=log(s);
                    med_ext= linterp(log_dist_AV,Aprior_VEC,lp,"linear");
                    std_ext = linterp(log_dist_AV,sigAprior_VEC,lp,"linear");
                    EF = GFunction(lav-med_ext,std_ext);
                    for(double mu=mag[0]-Mag_extN-nstd*err_mag[0];
                               mu<mag[0]-Mag_extN+nstd*err_mag[0];
                               mu+=dmu){
                        s = pow(10.,0.2*mu-2.);s2=s*s;lp=log(s);
                        // med_ext= linterp(log_dist_AV,Aprior_VEC,lp,"linear");
                        // std_ext = linterp(log_dist_AV,sigAprior_VEC,lp,"linear");
                        // EF = GFunction(lav-med_ext,std_ext);
                        tmp4=s2*s*tmp3*dmu*EF;
                        if(parallax_error>0.)
                            tmp4*=GFunction(parallax-1./s,parallax_error);
                        for(unsigned c=0;c<mag.size();++c){
                            Mag_ext=iso_grid->iso(i,j)->mag(k,mag_list[c])+av*mag_ext_const[c]*(1-av*avgradient[c]);
                            tmp4*=GFunction(mag[c]-(mu+Mag_ext),err_mag[c]);
                        }
                        if(tmp4==0.) continue;
                        if(prior->bprior()){
                            VecDoub X = conv::GalacticToCartesian({l,b,s});
                            tmp4*=prior->prior(X,
                                               iso_grid->fehgrid[i],
                                               agemodel);
                        }
                        total[0]+=tmp4;total[1]+=mu*tmp4;total[2]+=mu*mu*tmp4;
                        total[3]+=s*tmp4;total[4]+=s2*tmp4;
                        total[5]+=tmp4/s;total[6]+=tmp4/s2;
                        total[7]+=tmp4*log(agemodel);
                        total[8]+=tmp4*log(agemodel)*log(agemodel);
                        total[9]+=tmp4*mmodel;total[10]+=tmp4*mmodel*mmodel;
                        total[11]+=tmp4*iso_grid->fehgrid[i];total[12]+=tmp4*iso_grid->fehgrid[i]*iso_grid->fehgrid[i];
                        total[13]+=tmp4*lav;total[14]+=tmp4*lav*lav;
                        // Teff
                        total[15]+=tmp4*Teffmodel;
                        total[16]+=tmp4*Teffmodel*Teffmodel;
                        // logg
                        total[17]+=tmp4*loggmodel;
                        total[18]+=tmp4*loggmodel*loggmodel;
                        // Distane-Modulue-Age covariance
                        total[19]+=tmp4*log(agemodel)*mu;
                        if(return_pdf){
                            if(prob_mu.size()==0){
                                prob_mu.push_back(tmp4);
                                mu_grid.push_back(mu-grid_space*.5);
                                mu_grid.push_back(mu+grid_space*.5);
                            }
                            else{
                                if(mu<mu_grid.front()){
                                    while(mu<mu_grid.front()){
                                        mu_grid.insert(mu_grid.begin(),
                                                   mu_grid.front()-grid_space);
                                        prob_mu.insert(prob_mu.begin(),
                                                   tmp4);
                                    }
                                }
                                else if(mu>mu_grid.back()){
                                    while(mu>mu_grid.back()){
                                        mu_grid.push_back(mu_grid.back()+grid_space);
                                        prob_mu.push_back(tmp4);
                                    }
                                }
                                else{
                                    int bot,top;
                                    topbottom(mu_grid,mu,&bot,&top);
                                    prob_mu[bot]+=tmp4;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if(total[0]<=0.){
        if(NSTD>=20.){
            total[0]=std::numeric_limits<double>::infinity();
            total[1]=std::numeric_limits<double>::infinity();
            total[2]=-1.;
            for(int i=3;i<20;++i)
                total[i]=std::numeric_limits<double>::infinity();
            return total;
        }
        else return prob_distance_extinctprior(mag,Z,Teff,logg,l,b,err_mag,
                                               covar_ZTL,prior,return_pdf,
                                               2*NSTD,mag_list,Aprior_VEC,
                                               sigAprior_VEC,log_dist_AV,
                                               EL, parallax,parallax_error);
    }

    total=compute_mean_std_from_total(total);
    total[19]-=total[1]*total[7];
    if(return_pdf){
        double sum=0.;
        for(auto i:prob_mu)
            sum+=i*grid_space;
        for(unsigned i=0;i<prob_mu.size();++i)
            prob_mu[i]/=sum;
        prob_mu.push_back(mu_grid[0]);
        prob_mu.push_back(mu_grid[1]-mu_grid[0]);
        prob_mu.push_back(total[7]);
        prob_mu.push_back(total[8]);
        prob_mu.push_back(total[9]);
        prob_mu.push_back(total[10]);
        prob_mu.push_back(total[11]);
        prob_mu.push_back(total[12]);
        prob_mu.push_back(total[13]);
        prob_mu.push_back(total[14]);
        return prob_mu;
    }
    else
        return total;
}

template class DistanceCalculator<isochrone_johnson>;
template class DistanceCalculator<isochrone_padova>;
template class DistanceCalculator<isochrone_dartmouth>;
