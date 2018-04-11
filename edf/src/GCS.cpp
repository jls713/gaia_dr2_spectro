#include <Python.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "utils.h"
#include "GSLInterface/GSLInterface.h"
#include "coordtransforms.h"
#include "cubature/cubature.h"
#include "cuba/cuba.h"
#include "oct_int_exp.h"
using namespace std;
#include "tables_aa.h"
#include "potential.h"
#include "edf.h"
#include "sf.h"
#include "GR.h"
#include "GCS.h"

double dfsf(const VecDoub& v,double *sigma,double Vc,double F, double Fmean){
	return exp(-.5*(pow(v[0]/sigma[0],2)+pow((v[1]-Vc)/sigma[1],2)+pow(v[2]/sigma[2],2)+pow((F-Fmean)/sigma[3],2)));
}

#define NS 20

GCS_data::GCS_data(std::string IN){
	std::ifstream inFile;inFile.open(IN);
	if(!inFile.is_open()){std::cerr<<"Input file won't open."<<std::endl;}
	else std::cerr<<"Input file: "<<IN<<std::endl;
	std::string line; std::getline(inFile, line);
	VecDoub x_in (9,0), xerr(9,0);	std::string Name;
	while(inFile>>x_in[0]>>x_in[1]>>tmp>>x_in[2]>>x_in[3]>>x_in[4]
				>>x_in[5]>>x_in[6]>>x_in[7]>>x_in[8]
				>>xerr[3]>>xerr[4]
				>>xerr[5]>>xerr[6]>>xerr[7]>>xerr[8]){
			xerr[0]=0.;xerr[1]=0.;xerr[2]=0.12;
			InputCoords.push_back(x_in);InputCoords_e.push_back(xerr);
	}
	inFile.close();
	// Compute a fat grid of actions
	NDATA=InputCoords.size();
}

void GCS_data::find_actions(Potential_JS *Pot, Actions_AxisymmetricFudge_InterpTables *Tab, int readyForLL){
	std::ofstream outFile; outFile.open("GCSdata/GCSData_100414_withpolar.dat");
	outFile<<"R phi z vR vphi vz Z JR Lz Jz Rc"<<std::endl;
	VecDoub S;
	VecDoub SS = conv::StandardSolarPAUL;
	SS[3]-=239.1;SS[3]+=Pot->Vc(SS[0]);

	for(int i=0;i<NDATA;i++){
		VecDoub Gal = {InputCoords[i][3],InputCoords[i][4],1./InputCoords[i][5],InputCoords[i][6],InputCoords[i][7],InputCoords[i][8]};
		S = conv::GalacticToPolar(Gal,SS);
		VecDoub X ={S[0],0.,fabs(S[2]),S[3],-S[4],S[5]};
		if(Pot->H(X)>0) continue;//DF vanishes for E>0
		VecDoub Actions = Tab->actions(X);
		for(auto j: S)outFile<<j<<" ";
		outFile<<InputCoords[i][2]<<" ";
		// JR, Lz, Jz, Rc
		for(int i=0;i<4;i++)outFile<<Actions[i]<<" ";
		outFile<<std::endl;
	}
	outFile.close();

	if(readyForLL){
		VecDoub V(3,0); double Ferr,dist;
		ERROR_SAMPLES=50;
		NPOINTS=(int)(NDATA*ERROR_SAMPLES);
		for(int i=0; i<NDATA; i++){
			for(int j=0; j<ERROR_SAMPLES;j++){
				VecDoub X = {SS[0],0.,fabs(SS[1]),1e10,0.,0.};
				while(Pot->H(X)>0){
					// Scatter by errors
					for(int k=0;k<3;k++) V[k]=InputCoords[i][6+k]+rnGauss->nextnumber()*InputCoords_e[i][6+k];
					Ferr = InputCoords[i][2]+rnGauss->nextnumber()*InputCoords_e[i][2];
					// Find polar coords
					dist = -1.;
					while(dist<0.) dist = 1./(InputCoords[i][5]+InputCoords_e[i][5]*rnGauss->nextnumber());
					VecDoub Gal = {InputCoords[i][3],InputCoords[i][4],dist,V[0],V[1],V[2]};
					S = conv::GalacticToPolar(Gal,SS);
					X[0] = S[0]; X[2] = fabs(S[2]);
					X[3]=S[3];X[4]=-S[4];X[5]=S[5];
				}

				VecDoub Actions = Tab->actions(X);
				ActionTable.push_back(Actions);
				Metal.push_back(Ferr); Dist.push_back(dist);
			}
		}
		std::cout<<"Action tables computed"<<std::endl;
	}
	return;
}

void GCS_data::make_fake_table(sb15_edf_sf *edf, bool atSun, std::string outfile){
	VecDoub X; std::ofstream outFile; outFile.open(outfile);
	VecDoub max = find_max_values(edf, atSun);
	std::cerr<<max[0]<<" "<<max[1]<<" "<<max[2]<<std::endl;
	for(int i=0; i<NDATA; i++){
		X = Sampler(InputCoords[i][3], InputCoords[i][4], edf, InputCoords_e[i],max, atSun, outFile);
	}
	outFile.close();
}

double GCS_data::LogLike(sb15_edf_sf *edf){
	double Denom = log(edf->integrate_over_sphere(1e-3)), LogL=0.;
	// Error-free
	// for(int i=0; i<NDATA; i++){
	// 	VecDoub Gal = {InputCoords[i][3],InputCoords[i][4],1./InputCoords[i][5],InputCoords[i][6],InputCoords[i][7],InputCoords[i][8]};
	// 	VecDoub S = GalacticToPolar(Gal);
	// 	double x[2]={S[0],fabs(S[2])}; double v[3]={S[3]/100.,-S[4]/100.,S[5]/100.};
	// 	Phigl=Phi(x[0],x[1]); Phi0gl=Phi(x[0],0);
	// 	if(Phigl+.5*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2))>0) continue;
	// 	double dist = 1./InputCoords[i][5];
	// 	double LL = dist*dist*dist*dist*dist*dist*cos(InputCoords[i][4])*edf->chemDF_real_FULL(x,v,InputCoords[i][2],dist);
	// 	if(LL>0.) LogL+=log(LL);
	// }

	//With errors
//	edf->printParams();
	int n=0;
	for(int i=0; i<NDATA; i++){
		double sub = 0.;
		for(int k=0; k<ERROR_SAMPLES;k++){
			double d6 = Dist[n];d6*=d6;d6*=d6*d6;
			d6 = d6*cos(InputCoords[i][4])*edf->chemDF_actions(ActionTable[n], Metal[n], Dist[n]);
		 	sub += d6;
		 	n++;
		}
		if(sub>0.)LogL+=log(sub/(double)ERROR_SAMPLES);
	}
	return LogL-NDATA*Denom;
}


VecDoub find_max_values(sb15_edf_sf *edf, bool atSun){
	double f0;
	double Vc=edf->pot->Vc(edf->SunCoords()[0]);
	VecDoub X = {0.,0.,0.,0.,.7*Vc,0.};
	double FF, dist;
	double fmax=0,vmax=0.,Fmax=0.;
	double dv=.95*Vc/(double)(NS-1),dF=1./(double)(NS-1);
	double distmax = 0., dist_up = 0.45;
	if(atSun)dist_up=0.02;
	X[3]=0.;X[4]=.7*Vc;X[5]=0.;
	for(int j=0;j<NS;j++){//get peak value of DF
		FF=-0.5;
		for(int k=0;k<NS;k++){
			X[0]=edf->SunCoords()[0];X[2]=0.;
			for(dist=0.01;dist<dist_up;dist+=0.02){
				f0=dist*dist*edf->chemDF_real(X,FF,dist);
				if(f0>fmax){fmax=f0; vmax=X[4]; Fmax=FF; distmax=dist;}
			}
			FF+=dF;
		}
		X[4]+=dv;
	}
	VecDoub max = {vmax,Fmax,distmax};
	return max;
}

VecDoub Sampler(double l, double b, sb15_edf_sf *edf, VecDoub Err, VecDoub max, bool atSun, std::ofstream &outFile){
	// Samples a distance, metallicity and velocities at l,b
    double rr=1, rat=0;
	double f1,fs,sigma[4]={40.,40.,30.,0.4};
	double Vc=edf->pot->Vc(edf->SunCoords()[0]);
    VecDoub X ={edf->SunCoords()[0],0.,0.,0.,.7*Vc,0.};
	double FF, dist,distmax; VecDoub Gal(3), Pol(3);
	double fmax=0,vmax,Fmax;
	X[4]=max[0];vmax=max[0];Fmax=max[1];distmax=max[2];
	if(atSun) distmax=0.01;
	fmax = distmax*distmax*edf->chemDF_real(X,Fmax,distmax);

	dist = 0.01; // at Sun
	while(rr>rat){
		if(!atSun){
			dist = rn->nextnumber()*0.3; // sample a distance
			Gal = {l,b,dist};Pol = conv::GalacticToPolar(Gal,edf->SunCoords());
			X[0]=Pol[0];X[2]=fabs(Pol[2]);
		}
		else{
			Pol[0] = edf->SunCoords()[0]; Pol[1] = 0.00001; Pol[2]=edf->SunCoords()[1];
			X[0]=edf->SunCoords()[0];X[2]=edf->SunCoords()[1];
		}
		// Now sample a metallicity and a velocity
		for(int j=0;j<3;j++) X[j+3]=sigma[j]*rnGauss->nextnumber(); X[4]+=vmax;
		FF=sigma[3]*rnGauss->nextnumber()+Fmax;
		if(edf->pot->H(X)>0.) continue;//DF vanishes for E>0
		f1=dist*dist*edf->chemDF_real(X,FF,dist);
		fs=dfsf({X[3],X[4],X[5]},sigma,vmax,FF,Fmax);
		rat=f1/(1.8*fmax*fs);
		if(rat>1) std::cerr<<"RAT>1: "<<dist<<" "<<FF<<" "<<X[3]<<" "<<X[4]<<" "<<X[5]<<" "<<Fmax<<" "<<vmax<<" "<<rat<<" "<<f1<<" "<<fmax<<" "<<fs<<std::endl;
		rr=rn->nextnumber();
	}

	// find actions
	VecDoub Actions = edf->ActionCalculator->actions(X);
	double pthin = edf->thind->chemDF_actions(Actions, FF,dist);
	double pthick = edf->thickd->chemDF_actions(Actions, FF,dist);
	double phalo = edf->halod->chemDF_real(X, FF, dist);

	X[4]*=-1; // Important, as reverses direction of Galactic rotation
	VecDoub P = {Pol[0],Pol[1],Pol[2],X[3],X[4],X[5]};
	VecDoub G = conv::PolarToGalactic(P,edf->SunCoords());
	std::vector<double> VResults {G[2],G[3],G[4],G[5],FF};
	if(!atSun){
		double dtmp = -1.;
		while(dtmp<0.) dtmp = 1./(1./G[2]+rnGauss->nextnumber()*Err[5]);
		G[2] = dtmp;
		G[3] = (G[3]+rnGauss->nextnumber()*Err[6]);
		G[4] = (G[4]+rnGauss->nextnumber()*Err[7]);
		G[5] = (G[5]+rnGauss->nextnumber()*Err[8]);
	}
	FF = FF+rnGauss->nextnumber()*Err[2];
	P = conv::GalacticToPolar(G,edf->SunCoords());
	X[0] = P[0]; X[2] = fabs(P[2]);
	// find actions

	X[4]*=-1; // Important, as reverses direction of Galactic rotation
	Actions = edf->ActionCalculator->actions(X);
	for(auto i: P)outFile<<i<<" ";outFile<<FF<<" ";
	for(auto i: G)outFile<<i<<" ";
	for(int i=0;i<4;i++) outFile<<Actions[i]<<" ";
	outFile<<pthin<<" "<<pthick<<" "<<phalo<<std::endl;
	return VResults;
}
