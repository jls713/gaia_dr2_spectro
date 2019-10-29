#include "edf.h"
//===========================================================================*/
json edf_parameters(void){
    json parameters;
    std::string filename="config.json";
    std::ifstream inFile(filename);
    if(!inFile.good())
        std::cerr<<"Cannot open "<<filename<<std::endl;
    inFile>>parameters;
    inFile.close();
    return parameters;
}
//=============================================================================

double sb15_edf::full_DF(const VecDoub& X, double age, double RcP){
	if(pot->H(X)>0.) return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	return fullDF_actions(actions,age,RcP);
}

double sb15_edf::fullDF_actions(const VecDoub& actions, double age, double RcP){

	double Lz = actions[1],Jz=actions[2];
	double LcP=ActionCalculator->L_circ_current(RcP);
	double thin=0., thick=0.;
	double DV;
	// First thin
	if(age<thind->tau_T and RcP>0. and ifthin){
		DV = rm->DriftVelocity(age,RcP,Jz,thind->Rd(age));
		thin = GFunction(Lz-LcP-DV, rm->SigmaF(age, RcP,Jz))/rm->LzPrimeNormalization(LcP+DV, age, Jz);
		thin *= RcP*thind->fullDF_actions(actions, age, RcP);
	}
	if(age>thind->tau_T and RcP>0. and ifthick){
		DV = rm->DriftVelocity(age,RcP,Jz,
			             thickd->use_original_radius_for_drift?thickd->Rd(age):0.);
		thick = GFunction(Lz-LcP-DV, rm->SigmaF(age, RcP,Jz))/rm->LzPrimeNormalization(LcP+DV, age, Jz);
		thick *= RcP*thickd->fullDF_actions(actions, age, RcP);
	}
	double halo = 0.; // halo is meaningless when expressing in terms of birth radius
	return thin+thick+halo;
}

double sb15_edf::fullDF_actions_Z(const VecDoub& actions, double age, double Z){

	double Lz = actions[1], Jz=actions[2];
	double RcP = mr->RadiusFromMetal(age, Z);
	double thin=0., thick=0.;VecDoub ff;
	double DV;
	if(RcP>0. and Z>mr->FStart()){
		double LcP=ActionCalculator->L_circ_current(RcP);
		double dRdz = fabs(1./mr->dZdRc(RcP, age));
		// First thin
		if(age<thind->tau_T and ifthin){
			DV = rm->DriftVelocity(age,RcP,Jz,thind->Rd(age));
			thin= GFunction(Lz-LcP-DV, rm->SigmaF(age, RcP,Jz))/rm->LzPrimeNormalization(LcP+DV, 	age, Jz);
			thin *= RcP*dRdz*thind->fullDF_actions(actions, age, RcP);
		}
		if(age>thind->tau_T and ifthick){
			DV = rm->DriftVelocity(age,RcP,Jz,
			         thickd->use_original_radius_for_drift?thickd->Rd(age):0.);
			thick= GFunction(Lz-LcP-DV, rm->SigmaF(age, RcP,Jz))/rm->LzPrimeNormalization(LcP+DV, age, Jz);
			thick *= RcP*dRdz*thickd->fullDF_actions(actions, age, RcP);
		}
	}
	double halo=0.;
	if(age>thickd->tau_m-halod->delta_halo_age and ifhalo)
		halo = halod->chemDF_actions(actions,age,Z);
	return thin+thick+halo;
}

double sb15_edf::full_DF_Z(const VecDoub& X, double age, double Z){
	if(pot->H(X)>0.) return 0.;
	bool wf = true;
	VecDoub actions = ActionCalculator->actions(X,&wf);
	double Lz = actions[1], Jz=actions[2];
	double RcP = mr->RadiusFromMetal(age, Z);
	double thin=0., thick=0.;VecDoub ff;
	double DV;
	if(RcP>0. and Z>mr->FStart()){
		double LcP=ActionCalculator->L_circ_current(RcP);
		double dRdz = fabs(1./mr->dZdRc(RcP, age));
		// First thin
		if(age<thind->tau_T and ifthin){
			DV = rm->DriftVelocity(age,RcP,Jz,thind->Rd(age));
			thin= GFunction(Lz-LcP-DV, rm->SigmaF(age, RcP,Jz))/rm->LzPrimeNormalization(LcP+DV, 	age, Jz);
			thin *= RcP*dRdz*thind->fullDF_actions(actions, age, RcP);
		}
		if(age>thind->tau_T and ifthick){
			DV = rm->DriftVelocity(age,RcP,Jz,
			         thickd->use_original_radius_for_drift?thickd->Rd(age):0.);
			thick= GFunction(Lz-LcP-DV, rm->SigmaF(age, RcP,Jz))/rm->LzPrimeNormalization(LcP+DV, age, Jz);
			thick *= RcP*dRdz*thickd->fullDF_actions(actions, age, RcP);
		}
	}
	double halo=0.;
	if(age>thickd->tau_m-halod->delta_halo_age and ifhalo)
		halo = halod->chemDF_actions(actions,age,Z);
	return thin+thick+halo;
}
double sb15_edf::chemDF_actions(const VecDoub& Actions, double FeH){
	if(FeH<mr->FStart() and !ifhalo)return 0.; double sum=0.;
	if(ifthick)sum+=thickd->chemDF_actions(Actions,FeH);
	if(ifthin){sum+=thind->chemDF_actions(Actions,FeH);}
	if(ifhalo){sum+=halod->chemDF_actions(Actions,FeH);}
	return sum;
}

double sb15_edf::chemDF_real(const VecDoub& X, double FeH){
	if(FeH<mr->FStart() and !ifhalo)return 0.;
	if(pot->H(X)>0.) return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	double P = chemDF_actions(actions,FeH);
	if(P<0.)std::cerr<<"DF negative: Z="<<FeH<<std::endl;
	if(P!=P)std::cerr<<"DF nan: Z="<<FeH<<std::endl;
	return P;
}
double sb15_edf::chemDF_actions_justdiscs(const VecDoub& Actions, double FeH){
	if(FeH<mr->FStart())return 0.; double sum=0.;
	if(ifthick)sum+=thickd->chemDF_actions(Actions,FeH);
	if(ifthin){sum+=thind->chemDF_actions(Actions,FeH);}
	return sum;
}

double sb15_edf::chemDF_real_justdiscs(const VecDoub& X, double FeH){
	if(FeH<mr->FStart())return 0.;
	if(pot->H(X)>0.) return 0.;
	VecDoub actions = ActionCalculator->actions(X);
	double P = chemDF_actions_justdiscs(actions,FeH);
	if(P<0.){std::cerr<<"DF negative: Z="<<FeH<<" ";printVector(X);std::cout<<std::endl;}
	if(P!=P){std::cerr<<"DF nan: Z="<<FeH<<" ";printVector(X);std::cout<<std::endl;}
	return P;
}

//=============================================================================

sb15_edf::sb15_edf(Potential_JS *Pot, Actions_AxisymmetricFudge_InterpTables *Act, VecDoub Params, VecDoub StandardSolarInput, bool NewMetalRelation,bool ContinuousSFR,int order)
	:edf(Pot,Act)
	,rm(new RadialMigration(Params[9],tau_m,GAMMAT,Params[0],pot,StandardSolarInput[0],(Params.size()>15)?Params[15]:0.,(Params.size()>16)?Params[16]:0.))
	,mr(new tanh_metal_relation(Params[11],Params[13],Params[12],Params[10],tau_m))	,totalDiscMass(1.)
	,sfr(new SFR(T0,Params[8],tau_m, totalDiscMass))
	,StandardSolar(StandardSolarInput)
	,thind(new thin_disc_edf(Pot,Act,mr,sfr,rm,{Params[0],Params[1],Params[2],Params[3],tau_T,tau_1,BETA_R,BETA_Z},order))
	,thickd(new thick_disc_edf(Pot,Act,mr,sfr,rm,{Params[4],Params[5],Params[6],Params[7],tau_T,tau_m},order))
	,halod(new halo_edf(Pot,Act,{Params[14],FHALO,SIGFHALO}))
{
	StandardSolar[3]=Pot->Vc(StandardSolar[0])+StandardSolarInput[3];
	std::cout<<"Vc(Rsolar="<<StandardSolar[0]<<")="<<Pot->Vc(StandardSolar[0])<<std::endl;
	thind->setStandardSolar(StandardSolar);
	thickd->setStandardSolar(StandardSolar);
	rm->reset({Params[9],tau_m,GAMMAT,Params[0]},pot,StandardSolar[0],(Params.size()>15)?Params[15]:0.,(Params.size()>16)?Params[16]:0.);
	mr->reset({Params[11],Params[13],Params[12],Params[10],tau_m});
	sfr->reset({T0,Params[8],tau_m,totalDiscMass});
	thind->setParams({Params[0],Params[1],Params[2]*100.,Params[3]*100.,tau_T,tau_1,BETA_R,BETA_Z,0.,0.,0.,0.});
	thickd->setParams({Params[4],Params[5],Params[6]*100.,Params[7]*100.,tau_T,tau_m,0.,0.,0.,0.,0.});
	halod->setParams({Params[14],FHALO,SIGFHALO});
	ifthin=true;ifthick=true;ifhalo=false;
}

void sb15_edf::setParams(VecDoub Params,bool print){

	rm->reset({Params[9],tau_m,GAMMAT,Params[0]},pot,StandardSolar[0],(Params.size()>15)?Params[15]:0.,(Params.size()>16)?Params[16]:0.);
	mr->reset({Params[11],Params[13],Params[12],Params[10],tau_m});
	sfr->reset({T0,Params[8],tau_m,totalDiscMass});
	thind->setParams({Params[0],Params[1],Params[2],Params[3],tau_T,tau_1,BETA_R,BETA_Z,0.,0.,0.,0.});
	thickd->setParams({Params[4],Params[5],Params[6],Params[7],tau_T,tau_m,0.,0.,0.,0.,0.});
	halod->setParams({Params[14],FHALO,SIGFHALO,0.,0.,0.,0.,0.});
	if(print) printParams(std::cerr);
}

void sb15_edf::setParams(std::string jsonfile, bool print){
	// Much more complex setting of parameters from json file.
	// The use of a .edf file will set the 'default' parameters whilst
	// this interface allows for many more parameters to be set and various
	// options switched on and off.
    json p;
    std::ifstream inFile(jsonfile);
    if(!inFile.good())
        std::cerr<<"Cannot open "<<jsonfile<<std::endl;
    inFile>>p;
    inFile.close();
    rm->reset({p["rm"]["sigma_L"],p["sfr"]["tau_m"],GAMMAT,p["thin"]["Rd"]},pot,StandardSolar[0],
		p["rm"]["Jz0"],p["rm"]["sigR"]);

    // Change metal relation form to switch in json file
    delete mr;
    std::string metal_relation_type = p["metal_relation"];
    if(metal_relation_type=="tanh")
    	mr = new tanh_metal_relation(p["mr"]["fehstart"],p["mr"]["fehgrad"],p["mr"]["fehoffset"],p["mr"]["tauf"],p["sfr"]["tau_m"]);
    else if(metal_relation_type=="exp")
    	mr = new exp_metal_relation(p["mr"]["fehstart"],p["mr"]["fehgrad"],p["mr"]["fehoffset"],p["mr"]["tauf"],p["sfr"]["tau_m"]);
    else if(metal_relation_type=="linear")
    	mr = new linear_metal_relation(p["mr"]["fehstart"],p["mr"]["fehgrad"],p["mr"]["fehoffset"],p["mr"]["tauf"],p["sfr"]["tau_m"]);
    else
    	std::cout<<"Metal relation type not included in edf_params.json file\n";
    thind->set_mr(mr);thickd->set_mr(mr);
    // mr->reset({p["mr"]["fehstart"],p["mr"]["fehgrad"],p["mr"]["fehoffset"],p["mr"]["tauf"],p["sfr"]["tau_m"]});

	auto F = p["sfr"];
	if (F.find("type") != F.end()) {
		if(p["sfr"]["type"]=="interp"){
			delete sfr;
			VecDoub vals;
			for (auto it = p["sfr"].begin(); it != p["sfr"].end(); ++it)
				if(it.key()=="type" or it.key()=="tau_m" or it.key()=="tau_T")
					continue;
				else vals.push_back(it.value());
			int n=vals.size();
			sfr = new SFR_interp(vals,n,p["sfr"]["tau_m"]);
			thind->set_sfr(sfr);thickd->set_sfr(sfr);
		}
		else{
		    sfr->reset({p["sfr"]["T0"],p["sfr"]["tau_turn"],p["sfr"]["tau_m"],totalDiscMass});
		    T0 = p["sfr"]["T0"];
		}
	}
	else{
	    sfr->reset({p["sfr"]["T0"],p["sfr"]["tau_turn"],p["sfr"]["tau_m"],totalDiscMass});
	    T0 = p["sfr"]["T0"];
	}
    // Thin disc
    auto thin_json = p["thin"];
    VecDoub ThinParams={ thin_json["Rd"],	  thin_json["Rs"],
                      	 thin_json["sigmaR"], thin_json["sigmaZ"],
					  	 p["sfr"]["tau_T"]};
	auto extra_params = {"tau_1","BETA_R","BETA_Z","Rd_io","RSigmaZ","tau_1_z","tau_sigma"};
	VecDoub default_params = {tau_1,BETA_R,BETA_Z,0.,0.,0.};
	unsigned i=0;
	for(auto ex: extra_params){
		ThinParams.push_back(
		         (thin_json.find(ex) != thin_json.end())?(double)thin_json[ex]:default_params[i]);
		++i;
	}
	thind->setParams(ThinParams);
    // Thick disc
    auto thick_json = p["thick"];
    VecDoub ThickParams={thick_json["Rd"],	  thick_json["Rs"],
                      	 thick_json["sigmaR"],thick_json["sigmaZ"],
			 p["sfr"]["tau_T"],   p["sfr"]["tau_m"]};
	extra_params = {"RSigmaZ","Rd_io","sigmaR0","sigmaZ0","weight"};
	for(auto ex: extra_params)
		ThickParams.push_back(
		 (thick_json.find(ex) != thick_json.end())?(double)thick_json[ex]:0.);
	thickd->setParams(ThickParams);

	auto halo_json = p["halo"];
	VecDoub HaloParams = {halo_json["HW"],halo_json["FHALO"],halo_json["SIGFHALO"]};
	extra_params = {"J0","phi_coeff","z_coeff","halo_slope","delta_halo_age"};
	for(auto ex: extra_params)
		HaloParams.push_back(
		 (halo_json.find(ex) != halo_json.end())?(double)halo_json[ex]:0.);
	halod->setParams(HaloParams);
    if(thin_json.find("tau_1")!=thin_json.end())
        tau_1 = p["thin"]["tau_1"];
    tau_m = p["sfr"]["tau_m"];
    tau_T = p["sfr"]["tau_T"];
    BETA_R = p["thin"]["BETA_R"];
    BETA_Z = p["thin"]["BETA_Z"];
    SIGFHALO=p["halo"]["SIGFHALO"];
    FHALO=p["halo"]["FHALO"];
    if(p.find("use_average_radius") != p.end()){
    	if(p["use_average_radius"]==true){
	    	thind->set_average_radius(true);
	    	thickd->set_average_radius(true);
	    }
    }
    else{
    	thind->set_average_radius(false);
    	thickd->set_average_radius(false);
    }
    if(p.find("set_local_norm") != p.end()){
    	if(p["set_local_norm"]==true){
	    	thind->set_local_norm(true);
	    	thickd->set_local_norm(true);
	    }
    }
    else{
    	thind->set_local_norm(false);
    	thickd->set_local_norm(false);
    }
    if(p.find("use_original_radius") != p.end()){
    	if(p["use_original_radius"]==true){
	    	thickd->use_original_radius(true);
	    }
    }
    else{
    	thickd->use_original_radius(false);
    }
	if(print) printParams(std::cerr);
}

void sb15_edf::readParams(std::string ParamFile){
	std::ifstream inFile; inFile.open(ParamFile);VecDoub A(17,0);
	inFile>>A[0]>>A[1]>>A[2]>>A[3]>>A[4]>>A[5]>>A[6]>>A[7]>>A[8]>>A[9]
			>>A[10]>>A[11]>>A[12]>>A[13]>>A[14]>>A[15]>>A[16];
	A[9]*=100.;A[2]*=100.;A[3]*=100.;A[6]*=100.;A[7]*=100.;
	setParams(A,true);
	inFile.close();

}
void sb15_edf::printParams(std::ostream&o){
	o<<"Parameters:\nThin:\tR_d="<<thind->RD<<", R_s="<<thind->RSigma<<", s_R="<<thind->sigmaR<<", s_z="<<thind->sigmaZ;
	o<<"\nThick:\tR_d="<<thickd->RD<<", R_s="<<thickd->RSigma<<", s_R="<<thickd->sigmaR<<", s_z="<<thickd->sigmaZ<<", t_t="<<sfr->tau_turn;
	o<<"\nHalo:\t HW="<<halod->HW<<std::endl;
	o<<"Metallicity:\t FR="<<mr->FEHGrad<<", F0="<<mr->FEHOffset<<", F_m="<<mr->FEHStart<<", s_L="<<rm->SIGTHIN<<", gamma="<<rm->GAMMAT<<", tau_f="<<mr->TAUF<<", Jz0="<<rm->Jz0<<", sigR="<<rm->sigR<<std::endl<<std::endl;
}

std::string sb15_edf::params(void){
	std::string s("\nThin:\t");
	s+="R_d="+std::to_string(thind->RD)
	+", R_s="+std::to_string(thind->RSigma)
	+", s_R="+std::to_string(thind->sigmaR)
	+", s_z="+std::to_string(thind->sigmaZ);
	s+="\nThick:\t";
	s+="R_d="+std::to_string(thickd->RD)
	+", R_s="+std::to_string(thickd->RSigma)
	+", s_R="+std::to_string(thickd->sigmaR)
	+", s_z="+std::to_string(thickd->sigmaZ)
	+", t_t="+std::to_string(sfr->tau_turn);
	s+="\nHalo:\t";
	s+="HW="+std::to_string(halod->HW);
	s+="\nMetallicity:\t";
	s+="FR="+std::to_string(mr->FEHGrad)
	+", F0="+std::to_string(mr->FEHOffset)
	+", F_m="+std::to_string(mr->FEHStart)
	+", s_L="+std::to_string(rm->SIGTHIN)
	+", gamma="+std::to_string(rm->GAMMAT)
	+", tau_f="+std::to_string(mr->TAUF)
	+"\n\n";
	return s;
}

void sb15_edf::simpleprintParams(std::ofstream *output_file){
	VecDoub A = {	thind->RD,thind->RSigma,thind->sigmaR,thind->sigmaZ,
					thickd->RD,thickd->RSigma,thickd->sigmaR,thickd->sigmaZ,
					sfr->tau_turn, rm->SIGTHIN, mr->TAUF, mr->FEHStart, mr->FEHOffset, mr->FEHGrad, halod->HW};
	if(output_file) printVector_tofile(A, *output_file);
	printVector(A);
}
VecDoub sb15_edf::returnParams(void){
	VecDoub A = {thind->RD,thind->RSigma,thind->sigmaR,thind->sigmaZ,
					thickd->RD,thickd->RSigma,thickd->sigmaR,thickd->sigmaZ,
					sfr->tau_turn, rm->SIGTHIN, mr->TAUF, mr->FEHStart, mr->FEHOffset, mr->FEHGrad, halod->HW};
	return A;
}

//=============================================================================

const int nproc = 1;
const int nvec = 1; // maximum number of points per call
const int SEED = 0;//time(0);

double edf_integrate(integrand_t integrand, edf_norm_struct *P, double IE, double AE, std::string type, double *err, int verbosity, std::vector<VecDoub> *peaks){

    int neval,fail,nregions;
    double integral[1],error[1],prob[1];
    int NSIZE = P->x2min.size();
    double prod = 1.;
    for(int i=0;i<NSIZE;i++)prod*=(P->x2max[i]-P->x2min[i]);
    
    if(type=="Vegas")
        Vegas(NSIZE,nproc,integrand,P,nvec,IE,AE,verbosity,SEED,
        MINEVAL,MAXEVAL,NSTART,NINCREASE,NBATCH,GRIDNO,STATEFILE,SPIN,
        &neval,&fail,integral,error,prob);

    else if (type=="Suave")
        Suave(NSIZE,nproc,integrand,P,nvec,IE,AE,verbosity,SEED,
        MINEVAL,MAXEVAL,NNEW,FLATNESS,STATEFILE,SPIN,&nregions,
        &neval,&fail,integral,error,prob);

    else if (type=="Cuhre")
        Cuhre(NSIZE,nproc,integrand,P,nvec,IE,AE,verbosity,
        MINEVAL, MAXEVAL, 0, STATEFILE,SPIN,
        &nregions, &neval, &fail, integral, error, prob);

    else{
    	unsigned ngiven=0, ldxgiven=0;
        if(peaks){
            ngiven = peaks->size();
    	    ldxgiven = (*peaks)[0].size();
        }
    	double peaks_a[ldxgiven][ngiven];
    	for(unsigned l=0;l<ldxgiven;++l)
	    	for(unsigned n=0;n<ngiven;++n)
	    		peaks_a[l][n]=(*peaks)[n][l];
        Divonne(NSIZE,nproc,integrand,P,nvec,IE,AE,verbosity,SEED,
        MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
        BORDER, MAXCHISQ, MINDEVIATION,
        ngiven, ldxgiven, *peaks_a, NEXTRA, nullptr,STATEFILE,SPIN,
        &nregions, &neval, &fail, integral, error, prob);
    }
    if(err)*err=prod*error[0];
    if(fail!=0)std::cerr<<"Error: Required accuracy not reached:" <<fail<<std::endl;
    return prod*integral[0];
}

double edf_integrate_peak(integrand_t integrand, edf_norm_struct *P, VecDoub peak, double IE, double AE, double *err){

    int neval,fail,nregions;
    double integral[1],error[1],prob[1];
    int NSIZE = P->x2min.size();
    double prod = 1.;
    for(int i=0;i<NSIZE;i++)prod*=(P->x2max[i]-P->x2min[i]);
    int ngiven=1, ldxgiven=peak.size();

    Divonne(NSIZE,nproc,integrand,P,1,IE,AE,0,SEED,
        MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
        BORDER, MAXCHISQ, MINDEVIATION,
        ngiven, ldxgiven, &peak[0], NEXTRA, nullptr,STATEFILE,SPIN,
        &nregions, &neval, &fail, integral, error, prob);
    if(err)*err=prod*error[0];
    if(fail!=0)std::cerr<<"Error: Required accuracy not reached:" <<fail<<std::endl;
    return prod*integral[0];
}
static int numericalnorm_integrand_cuba(const int ndim[],const double s[], const int*fdim, double fval[], void *fdata){
    VecDoub actions(3,0.); double prod=1.;
    edf_norm_struct *P = (edf_norm_struct *) fdata;
    for(int i=0;i<3;i++){
    	if(s[i]==1. or s[i]==0.){fval[0]=0.; return 0;}
    	prod*=(1.-s[i])*(1.-s[i]);
    	actions[i]=P->scale*s[i]/(1.-s[i]);
    }
    actions.push_back(P->edf->ActionCalculator->R_circ(actions[1]));
    VecDoub f = P->edf->ActionCalculator->PotentialFreq(actions[3]);
    for(auto i:f) actions.push_back(i);
    fval[0]=P->edf->chemDF_actions(actions,s[3])/prod;
	if(fval[0]!=fval[0]) fval[0]=0.;
    return 0;
}

double sb15_edf::normalization(double scale){
	VecDoub x2min(4,0.), x2max(4,1.);
	x2max[3]=mr->FeH(0.,0.)-1e-5;x2min[3]=mr->FStart()+1e-5;
	edf_norm_struct P(this,x2min,x2max,scale);
	return pow(2.*PI*scale,3)*edf_integrate(&numericalnorm_integrand_cuba,&P,1e-4,0,"Divonne");
}
static int density_tt_integrand_cuba(const int ndim[], const double y[], const int *fdim, double fval[], void *fdata){
	double y2[4];
	edf_density_struct *P = (edf_density_struct *) fdata;
	for(int i=0;i<4;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
	VecDoub X = {P->R,0.,P->Z,y2[0], y2[1], y2[2]};
	if(P->edf->pot->H(X)>P->edf->pot->Phi({1.e7*(X[0]==0.?1.:X[0]),1.e7*(X[1]==0.?1.:X[1]),1.e7*(X[2]==0.?1.:X[2])}))
		{ fval[0]=0.;return 0;}
        fval[0] = P->edf->chemDF_real_justdiscs(X,y2[3]);
	return 0;
}

static int density_halo_integrand_cuba(const int ndim[], const double y[], const int *fdim, double fval[], void *fdata){
	double y2[4];
	edf_density_struct *P = (edf_density_struct *) fdata;
	for(int i=0;i<4;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
	VecDoub X = {P->R,0.,P->Z,y2[0], y2[1], y2[2]};
	if(P->edf->pot->H(X)>P->edf->pot->Phi({1.e7*(X[0]==0.?1.:X[0]),1.e7*(X[1]==0.?1.:X[1]),1.e7*(X[2]==0.?1.:X[2])}))
		{ fval[0]=0.;return 0;}
	fval[0] = P->edf->halod->chemDF_real(X,y2[3]);
	return 0;
}

double sb15_edf::density(double R, double z, double IE, VecDoub Flims){

	double Ve=sqrt(-2*(pot->Phi({R,0.,z})-pot->Phi({10*R,0.,10*z}))), Vze=Ve;
	double Flim=mr->FeH(0.,0.)-1e-5;
	VecDoub Flimits = Flims;
	if(Flims[0]<mr->FStart()) Flimits[0]=mr->FStart()+1e-5;
	if(Flims[1]>Flim) Flimits[1]=Flim;

	VecDoub x2min={-.8*Ve,0.01*Ve,0.,Flimits[0]};
	VecDoub x2max={.8*Ve,.8*Ve,.8*Vze,Flimits[1]};

	edf_density_struct P(this,R,z,x2min,x2max);

  	double thinthick=0.,halobit=0.;
        if(ifthin or ifthick)
		thinthick=2.*edf_integrate(&density_tt_integrand_cuba,&P,IE,0);

	if(!ifhalo) return thinthick;

	if(Flims[0]<-3.5) Flimits[0]=-3.5;
	if(Flims[1]>.5) Flimits[1]=0.5;

	VecDoub x2min2={-.8*Ve,-.8*Ve,-.8*Vze,Flimits[0]};
	VecDoub x2max2={.8*Ve,.8*Ve,.8*Vze,Flimits[1]};

	edf_density_struct P2(this,R,z,x2min2,x2max2);

	halobit=edf_integrate(&density_halo_integrand_cuba,&P2,IE,0);
	
        return thinthick+halobit;
}

static int FHistIntegrand_cuba(const int ndim[], const double y[], const int *fdim, double fval[], void *fdata){
	double y2[3];
	edf_vhist_struct *P = (edf_vhist_struct *) fdata;
	for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
	VecDoub X = {P->R,0.,P->Z,y2[0], y2[1], y2[2]};
	fval[0] = P->edf->chemDF_real_justdiscs(X,P->V)+2.*P->edf->halod->haloDF_nometal_real(X);
	return 0;
	}

static int UHistIntegrand_cuba(const int ndim[], const double y[], const int *fdim, double fval[], void *fdata){
	double y2[3];
	edf_vhist_struct *P = (edf_vhist_struct *) fdata;
	for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
	VecDoub X = {P->R,0.,P->Z,P->V, y2[0], y2[1]};
	fval[0] = P->edf->chemDF_real_justdiscs(X,y2[2])+2.*P->edf->halod->haloDF_nometal_real(X)/(P->x2max[2]-P->x2min[2]);
	return 0;
	}

static int VHistIntegrand_cuba(const int ndim[], const double y[], const int *fdim, double fval[], void *fdata){
	double y2[3];
	edf_vhist_struct *P = (edf_vhist_struct *) fdata;
	for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
	VecDoub X = {P->R,0.,P->Z, y2[0],P->V, y2[1]};
	fval[0] = P->edf->chemDF_real_justdiscs(X,y2[2])+2.*P->edf->halod->haloDF_nometal_real(X)/(P->x2max[2]-P->x2min[2]);
	return 0;
	}

static int WHistIntegrand_cuba(const int ndim[], const double y[], const int *fdim, double fval[], void *fdata){
	double y2[3];
	edf_vhist_struct *P = (edf_vhist_struct *) fdata;
	for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
	VecDoub X = {P->R,0.,P->Z, y2[0], y2[1],P->V};
	fval[0] = P->edf->chemDF_real_justdiscs(X,y2[2])+2.*P->edf->halod->haloDF_nometal_real(X)/(P->x2max[2]-P->x2min[2]);
	return 0;
	}

struct edf_pmhist_struct: edf_norm_struct{
	double l,b,s,r,phi,z;
	int which; // which=0 -> pm_l, which=1 -> pm_b
	edf_pmhist_struct(sb15_edf *edf, double l, double b, double s, VecDoub x2m, VecDoub x2n,int which)
		:edf_norm_struct(edf,x2m,x2n,1.),l(l),b(b),s(s),which(which){
			VecDoub P = conv::GalacticToPolar({l,b,s},edf->SunCoords());
			r=P[0];phi=P[1];z=P[2];
		}
};

static int PMIntegrand_cuba(const int ndim[], const double y[], const int *fdim, double fval[], void *fdata){
	double y2[4];
	edf_pmhist_struct *P = (edf_pmhist_struct *) fdata;
	for(int i=0;i<4;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
	VecDoub X = {P->r,P->phi,P->z, y2[0], -y2[1], y2[2]};
	// if(P->edf->pot->H(X)>P->edf->pot->Phi({1.e7*(X[0]==0.?1.:X[0]),1.e7*(X[1]==0.?1.:X[1]),1.e7*(X[2]==0.?1.:X[2])}))
	// 	{ fval[0]=0.;return 0;}
	VecDoub G = conv::PolarToGalactic(X,P->edf->SunCoords());
	double ff = G[4]; if(P->which==1) ff = G[5];
	X = {P->r,0.,P->z, y2[0], y2[1], y2[2]};
	fval[0] = ff*P->edf->chemDF_real_justdiscs(X,y2[3]);
	return 0;
	}

static int PMIntegrand_Halo_cuba(const int ndim[], const double y[], const int *fdim, double fval[], void *fdata){
	double y2[4];
	edf_pmhist_struct *P = (edf_pmhist_struct *) fdata;
	for(int i=0;i<4;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
	VecDoub X = {P->r,P->phi,P->z, y2[0], -y2[1], y2[2]};
	if(P->edf->pot->H(X)>P->edf->pot->Phi({1.e7*(X[0]==0.?1.:X[0]),1.e7*(X[1]==0.?1.:X[1]),1.e7*(X[2]==0.?1.:X[2])}))
		{ fval[0]=0.;return 0;}
	VecDoub G = conv::PolarToGalactic(X,P->edf->SunCoords());
	double ff = G[4]; if(P->which==1) ff = G[5];
	X = {P->r,0.,P->z, y2[0], y2[1], y2[2]};
	fval[0] = ff*P->edf->halod->chemDF_real(X,y2[3]);
	return 0;
	}

double sb15_edf::FHist(double R, double z, double F, double IE){

	double Ve=sqrt(-2*(pot->Phi({R,0.,z})-pot->Phi({10*R,0.,10*z}))), Vze=Ve;
	VecDoub x2min={0.01*Ve,0.01*Ve,-.8*Vze};
	VecDoub x2max={0.8*Ve,.8*Ve,.8*Vze};

	edf_vhist_struct P(this,R,z,F,x2min,x2max);

	int neval, fail;
  	double integral[1], error[1], prob[1];

	Vegas(3,1,FHistIntegrand_cuba,&P,1,IE,EPSABS,0,0,
	MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE,SPIN,
	&neval, &fail, integral, error, prob);

    double prod=1;
    for(int i=0; i<3; i++)prod*=(x2max[i]-x2min[i]);
	return 2.*prod*integral[0];
}

double sb15_edf::UHist(double R, double z, double v, double IE){

	double Ve=sqrt(-2*(pot->Phi({R,0.,z})-pot->Phi({10*R,0.,10*z}))), Vze=Ve;
	double Flim=mr->FeH(0.,0.)-0.00001;
	VecDoub x2min={0.001*Ve,-.8*Vze,mr->FStart()+1e-5};
	VecDoub x2max={.8*Ve,.8*Vze,Flim};
	edf_vhist_struct P(this,R,z,v,x2min,x2max);
	return edf_integrate(&UHistIntegrand_cuba,&P,IE);
}

double sb15_edf::VHist(double R, double z, double v, double IE){

	double Ve=sqrt(-2*(pot->Phi({R,0.,z})-pot->Phi({10*R,0.,10*z}))), Vze=Ve;
	double Flim=mr->FeH(0.,0.)-0.00001;
	VecDoub x2min={0.,-.8*Vze,mr->FStart()+1e-5};
	VecDoub x2max={.8*Ve,.8*Vze,Flim};
	edf_vhist_struct P(this,R,z,v,x2min,x2max);
	return edf_integrate(&VHistIntegrand_cuba,&P,IE);
}

double sb15_edf::WHist(double R, double z, double v, double IE){

	double Ve=sqrt(-2*(pot->Phi({R,0.,z})-pot->Phi({10*R,0.,10*z})));
	double Flim=mr->FeH(0.,0.)-mr->FStart();
	VecDoub x2min={-.8*Ve,0.01*Ve,mr->FStart()+1e-5};
	VecDoub x2max={.8*Ve,.8*Ve,Flim};
	edf_vhist_struct P(this,R,z,v,x2min,x2max);
	return edf_integrate(&WHistIntegrand_cuba,&P,IE);
}

static int LOSIntegrand_cuba(const int ndim[], const double y[], const int *fdim, double fval[], void *fdata){
	double y2[5];
	edf_density_struct *P = (edf_density_struct *) fdata;
	for(int i=0;i<5;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
	VecDoub Gal = {P->R,P->Z,y2[0]};
	VecDoub Pol = conv::GalacticToPolar(Gal,P->edf->SunCoords());
	VecDoub X = {Pol[0],0.,fabs(Pol[2]),y2[1],y2[2],y2[3]};
	fval[0] = P->edf->chemDF_real(X,y2[4])+4.*P->edf->halod->haloDF_nometal_real(X)/(P->x2max[4]-P->x2min[4]);
	return 0;
	}

double sb15_edf::integrate_along_line_of_sight(double l, double b, double IE){
	// Doesn't include halo contribution properly
	double R = StandardSolar[0];double z=0.;
	double Ve=sqrt(-2*(pot->Phi({R,0.,z})-pot->Phi({10*R,0.,10*z}))), Vze=Ve;
	double Flim=mr->FeH(0.,0.)-0.00001;
	VecDoub x2min={0.,-.8*Ve,0.01*Ve,0.,mr->FStart()};
	VecDoub x2max={0.5,.8*Ve,.8*Ve,.8*Vze,Flim};
	edf_density_struct P(this,l,b,x2min,x2max);
	return edf_integrate(&LOSIntegrand_cuba,&P,IE);
}

VecDoub sb15_edf::mean_proper_motions(double l, double b, double s, double IE){
	VecDoub PP = conv::GalacticToPolar({l,b,s},SunCoords());
	double Ve=sqrt(-2*(pot->Phi(PP)-pot->Phi(PP*10.))), Vze=Ve;
	double Flim=mr->FeH(0.,0.)-0.00001;
	VecDoub x2min={-.8*Ve,0.01*Ve,-.8*Vze,mr->FStart()};
	VecDoub x2max={.8*Ve,.8*Ve,.8*Vze,Flim};
	edf_pmhist_struct PM_struct(this,l,b,s,x2min,x2max,0);
	double fpml = edf_integrate(&PMIntegrand_cuba,&PM_struct,IE);
	edf_pmhist_struct PM_struct2(this,l,b,s,x2min,x2max,1);
	double fpmb = edf_integrate(&PMIntegrand_cuba,&PM_struct2,IE);
	edf_density_struct Dens_struct(this,PP[0],PP[2],x2min,x2max);
	double dens = 2.*edf_integrate(&density_tt_integrand_cuba,&Dens_struct,IE,0);
	// VecDoub Flimits(2,0.);
	// Flimits[0]=-3.5;
	// Flimits[1]=0.5;

	// VecDoub x2min2={-.8*Ve,-.8*Ve,-.8*Vze,Flimits[0]};
	// VecDoub x2max2={.8*Ve,.8*Ve,.8*Vze,Flimits[1]};
	// edf_pmhist_struct PM_Halo_struct(this,l,b,s,x2min2,x2max2,0);
	// fpml += edf_integrate(&PMIntegrand_Halo_cuba,&PM_Halo_struct,IE);
	// edf_pmhist_struct PM_Halo_struct2(this,l,b,s,x2min2,x2max2,0);
	// fpmb += edf_integrate(&PMIntegrand_Halo_cuba,&PM_Halo_struct2,IE);
	// edf_density_struct Dens_Halo_struct(this,PP[0],PP[2],x2min2,x2max2);
	// dens += edf_integrate(&density_halo_integrand_cuba,&Dens_Halo_struct,IE,0);

	return {fpml/dens,fpmb/dens};
}

static int moments_Integrand_cuba(const int ndim[], const double y[], const int *fdim, double fval[], void *fdata){
	double y2[3];
	edf_moments_struct *P = (edf_moments_struct *) fdata;
	for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];
	VecDoub X = {P->R,0.,P->Z,y2[0], y2[1], y2[2]};
	fval[0] = P->edf->full_DF_Z(X,P->tau,P->met);
	fval[1] = y2[1]*fval[0];// mean Vphi
	for(int i=0;i<3;++i)
		fval[i+2] = y2[i]*y2[i]*fval[0];// disp
	return 0;
	}

VecDoub edf_multi_integrate(integrand_t integrand, int NINTS, edf_norm_struct *P, double IE, double AE, std::string type, double *err, int verbosity){

    int neval,fail,nregions;
    double integral[NINTS],error[NINTS],prob[NINTS];
    int NSIZE = P->x2min.size();
    double prod = 1.;
    for(int i=0;i<NSIZE;i++)prod*=(P->x2max[i]-P->x2min[i]);

    if(type=="Vegas")
        Vegas(NSIZE,NINTS,integrand,P,nvec,IE,AE,verbosity,SEED,
        MINEVAL,MAXEVAL,NSTART,NINCREASE,NBATCH,GRIDNO,STATEFILE,SPIN,
        &neval,&fail,integral,error,prob);

    else if (type=="Suave")
        Suave(NSIZE,NINTS,integrand,P,nvec,IE,AE,verbosity,SEED,
        MINEVAL,MAXEVAL,NNEW,FLATNESS,STATEFILE,SPIN,&nregions,
        &neval,&fail,integral,error,prob);

    else if (type=="Cuhre")
        Cuhre(NSIZE,NINTS,integrand,P,nvec,IE,AE,verbosity,
        MINEVAL, MAXEVAL, 0, STATEFILE,SPIN,
        &nregions, &neval, &fail, integral, error, prob);

    else
        Divonne(NSIZE,NINTS,integrand,P,nvec,IE,AE,verbosity,SEED,
        MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
        BORDER, MAXCHISQ, MINDEVIATION,
        NGIVEN, LDXGIVEN, nullptr, NEXTRA, nullptr,STATEFILE,SPIN,
        &nregions, &neval, &fail, integral, error, prob);
    if(err)*err=prod*error[0];
    if(fail!=0)std::cerr<<"Error: Required accuracy not reached:" <<fail<<std::endl;
    VecDoub results(NINTS,0.);
    for(unsigned i=0;i<NINTS;++i)
    	results[i]=prod*integral[i];
    return results;
}

VecDoub sb15_edf::moments(double R, double z, double tau, double Z, double IE){
	double Ve=sqrt(-2*(pot->Phi({R,0.,z})-pot->Phi({10*R,0.,10*z}))), Vze=Ve;
	VecDoub x2min={-.9*Ve,0.01*Ve,-.9*Vze};
	VecDoub x2max={.9*Ve,.9*Ve,.9*Vze};
	edf_moments_struct P(this,R,z,Z,tau,x2min,x2max);
	VecDoub results = edf_multi_integrate(&moments_Integrand_cuba,5,&P,IE);
	for(unsigned i=1;i<5;++i)
		results[i]/=results[0];
	results[3]=sqrt(results[3]-results[1]*results[1]);
	results[2]=sqrt(results[2]);
	results[4]=sqrt(results[4]);
	return results;
}
