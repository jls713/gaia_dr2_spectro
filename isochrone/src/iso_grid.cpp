#include "iso_grid.h"
//=============================================================================

VecDoub load_padova_metalorage_file(std::string filename){
    std::ifstream infile;infile.open(filename);
    VecDoub Zout; double Z;
    while(infile>>Z)
        Zout.push_back(Z);
    infile.close();
    return Zout;
}

template<class isochrone_g>
void isochrone_grid<isochrone_g>::thin_feh_grid(double feh_err){
    double mean_spacing=(fehgrid.back()-fehgrid.front())/(double)(fehgrid.size());
    while(mean_spacing<feh_err/2.){
        for(unsigned i=1;i<fehgrid.size()-(fehgrid.size()%2==1);++i)
            fehgrid.erase(fehgrid.begin()+i);
        mean_spacing=(fehgrid.back()-fehgrid.front())/(double)(fehgrid.size());
    }
}

//=============================================================================

template<class isochrone_g>
isochrone_grid<isochrone_g>::isochrone_grid(std::string type,int thin, double feh_err)
    {
    // find list of isochrone files
    std::vector<std::string> iso_files;
    std::string iso_dir = parameters()["dir"]["isochrones"];
    std::string dir;

    if(type=="BaSTI"){
        dir = iso_dir+"BaSTI_Johnson/more";
        agegrid = BaSTI_More::ageList;
        fehgrid = BaSTI_More::FeHList;
    }
    else if(type=="Padova"){
        dir = iso_dir+"PARSEC_Gaia/grid/";
        agegrid = load_padova_metalorage_file(iso_dir+
                                              "PARSEC_Gaia/age_vals.dat");
        fehgrid = load_padova_metalorage_file(iso_dir+
                                              "PARSEC_Gaia/metal_vals.dat");
    }
    else if(type=="Dartmouth"){
        dir = iso_dir+"Dartmouth/UBVRIJHKsKp/afe0";
        agegrid = Dartmouth::ageList;fehgrid = Dartmouth::FeHList;
    }
    NF = fehgrid.size(); NA = agegrid.size();
    iso_grid_size = NA*NF;
    GetFilesInDirectory(iso_files, dir);
    int isoPadova_filters=5;
    if(iso_files.size()!=iso_grid_size and (type=="BaSTI"))
        std::cout<<"More isochrones ("<<iso_files.size()<<") than grid spaces("
            <<iso_grid_size<<")!"<<std::endl;
    if(iso_files.size()!=iso_grid_size/NA*isoPadova_filters and (type=="Padova"))
        std::cout<<"More isochrones ("<<iso_files.size()/isoPadova_filters<<") than grid spaces("
            <<iso_grid_size/NA<<")!"<<std::endl;
    if(iso_files.size()!=iso_grid_size/NA*2 and (type=="Dartmouth"))
        std::cout<<"More isochrones ("<<iso_files.size()<<") than grid spaces("
            <<iso_grid_size/NA*2<<")!"<<std::endl;

    int iso_grid_size_orig = iso_grid_size;
    thin_feh_grid(feh_err);
    NF = fehgrid.size(); NA = agegrid.size();
    iso_grid_size = NA*NF;

    iso_grid = std::vector<isochrone_g>(iso_grid_size);

    for(unsigned i=0;i<NF;++i)
        massmax.push_back(VecDoub(NA,0.));

    // fill iso_grid with isochrones
    int i, j; double f, a;
    if(type=="Dartmouth"){
        auto filters = parameters()["isochrones"]["filters"];
        bool johnson=false;
        for(auto f: filters){
            if(f!="Johnson")
                std::cerr<<f<<" filter set not available for Dartmouth."<<std::endl;
            else
                johnson=true;
        }
        if(johnson==false)
            std::cerr<<"Johnson filter set not set in config.json but is the only filter set available for Dartmouth so loading anyway\n";
        for(auto file: iso_files){
            VecDoub AA = Dartmouth::ageList_2;
            if(file.back()=='2') AA = Dartmouth::ageList_1;
            for(auto age_index:AA){
                isochrone_g iso;
                iso.fill({file}, dir, age_index);
                // sort results
                for(i=0;i<NF;i++){
                    f = fehgrid[i];if(fabs(f-iso.feh())<0.0001) break;}
                for(j=0;j<NA;j++){
                    a = agegrid[j];if(fabs(a-iso.tau())<0.0001) break;}
                if(fabs(f-iso.feh())<0.0001){
                    iso_grid[j+NA*i]=iso;
                    massmax[i][j]=iso.maxmass();
                }
            }
        }
    }
    else if(type=="Padova"){
	    int i, j;
	    double f, a;
        auto filters = parameters()["isochrones"]["filters"];
	    isoPadova_filters=filters.size();
        std::vector<std::vector<std::string>> isoch(iso_grid_size_orig/NA,std::vector<std::string>(isoPadova_filters,""));
	    int ny=0;
	    std::map<std::string,std::string> filter_codes={{"Gaia","g"},
                                                        {"SDSS","s"},
                                                        {"Johnson","u"},
							                            {"2MASS","2"},
                                                        {"PanSTARRS","p"}};
	    for(std::string str:{"g","u","2","s","p"}){
            bool cont=true;
    		for(auto f: filters)
                if (str.compare(0,1,filter_codes[f])==0
                    or filters[0]=="All") cont=false;
            if(cont)continue;
            auto it = std::find_if(iso_files.begin(),iso_files.end(),
                                   [str,dir](const std::string& obj){
				return obj.compare(dir.length()+1,1,str)==0;});
    		int nn=0;
    		while(it!=std::end(iso_files)){
    			isoch[nn][ny]=*it;
    			++nn;
    			it = std::find_if(std::next(it),iso_files.end(),
                                  [str,dir](const std::string& obj){
    				return obj.compare(dir.length()+1,1,str)==0;});
    		}
    		++ny;
	    }
	    for(auto s: isoch){
            isochrone_g iso_t;
            double ff = iso_t.get_metallicity(s,dir);
            for(i=0;i<NF;i++){
                f = fehgrid[i];
                if(fabs(f-ff)<0.0001) break;
            }
            if(fabs(f-ff)>0.0001)
                continue;
            for(auto age_index:agegrid){
    		    isochrone_g iso;
    		    iso.fill(s, dir, age_index);
    		    // sort results
    		    for(i=0;i<NF;i++){
        			f = fehgrid[i];
                    if(fabs(f-iso.feh())<0.0001) break;
                }
	            for(j=0;j<NA;j++){
        			a = agegrid[j];
                    if(fabs(a-iso.tau())<0.0001) break;
                }
	            if(fabs(f-iso.feh())<0.0001){
    			    iso_grid[j+NA*i]=iso;
    			    massmax[i][j]=iso.maxmass();
    		    }
    		}
            std::cout<<iso_grid[j+NA*i].feh()<<std::endl;
	    }
   }
    else{
        bool johnson=false;
        auto filters = parameters()["isochrones"]["filters"];
        for(auto f: filters){
            if(f!="Johnson"){
                std::cerr<<f<<" filter set not available for BaSTI. ";
                std::cerr<<"If Sloan, this is possible but not implemented.\n";
            }
            else
                johnson=true;
        }
        if(johnson==false){
            std::cerr<<"Johnson filter set not set in config.json but is ";
            std::cerr<<"the only filter set available for BaSTI so loading anyway\n";
        }
        for(auto file: iso_files){
            isochrone_g iso;
            double ff = iso.get_metallicity({file},dir);
            for(i=0;i<NF;i++){
                f = fehgrid[i];
                if(fabs(f-ff)<0.0001) break;
            }
            if(fabs(f-ff)>0.0001)
                continue;
            iso.fill({file}, dir, 0.);
            // sort results
            for(i=0;i<NF;i++){
                f = fehgrid[i];if(fabs(f-iso.feh())<0.0001) break;}
            for(j=0;j<NA;j++){
                a = agegrid[j];if(fabs(a-iso.tau())<0.0001) break;}
            if(fabs(f-iso.feh())<0.0001){
                    iso_grid[j+NA*i]=iso;
                    massmax[i][j]=iso.maxmass();
            }
        }
    }
    for(unsigned i=0;i<NF;++i)
    for(unsigned j=0;j<NA;++j)
        for(unsigned m=0;m<iso(i,j)->N();++m)
            iso(i,j)->set_maxage(m,max_age(iso(i,j)->feh(),iso(i,j)->initial_mass(m)));
    std::cerr<<type+" Isochrones loaded"<<std::endl;
    auto filters = parameters()["isochrones"]["filters"];
    std::cerr<<"Using filters: ";
        for(auto f: filters)
            std::cerr<<f<<", ";
    std::cerr<<std::endl;
    std::cerr<<"Min age: "<<agegrid[0]<<", Max age: "<<agegrid[NA-1]<<", "<<NA<<" ages."<<std::endl;
    std::cerr<<"Min feh: "<<fehgrid[0]<<", Max feh: "<<fehgrid[NF-1]<<", "<<NF<<" metallicities."<<std::endl;
}
//=============================================================================

template<class isochrone_g>
isochrone_g* isochrone_grid<isochrone_g>::iso(int i,int j){
    if(j>NA-1 or i>NF-1) std::cerr<<"Outside isochrone grid"<<std::endl;
    return &iso_grid[j+NA*i];
}
//=============================================================================

template<class isochrone_g>
std::vector<int> isochrone_grid<isochrone_g>::find_nearest(double Z, double age, double M){

    std::vector<int> V(3,0);
    int bot_Z,top_Z,bot_a,top_a,bot_m=0,top_m=1;

    if(Z<fehgrid.front()) V[0]=0;
    else if(Z>fehgrid.back()) V[0]=fehgrid.size()-1;
    else{
    	topbottom(fehgrid,Z,&bot_Z,&top_Z,"Isochrone Metal");
    	if(Z-fehgrid[bot_Z]>fehgrid[top_Z]-Z)
    	    V[0]=top_Z;
    	else
    	    V[0]=bot_Z;
    }
	if(age<agegrid.front()) V[1]=0;
    else if(age>agegrid.back()) V[1]=agegrid.size()-1;
    else{
    	topbottom(agegrid,age,&bot_a,&top_a,"Isochrone Age");
    	if(age-agegrid[bot_a]>agegrid[top_a]-age)
    	    V[1]=top_a;
    	else
    	    V[1]=bot_a;
    }

    int m;
    top_m=iso(V[0],V[1])->N()-1;
    if(M<iso(V[0],V[1])->initial_mass(bot_m)) V[2]=-1;
    else if(M>iso(V[0],V[1])->initial_mass(top_m)) V[2]=-1; // Unphysical
    else{
    	while(top_m-bot_m>1){
    	    m=(top_m+bot_m)/2;
    	    if((iso(V[0],V[1])->initial_mass(top_m)-M)*(M-iso(V[0],V[1])->initial_mass(m))>=0.) 	bot_m=m;
    	    else top_m=m;
    	}
    	if(M-iso(V[0],V[1])->initial_mass(bot_m)>iso(V[0],V[1])->initial_mass(top_m)-M)
    	    V[2]=top_m;
    	else
    	    V[2]=bot_m;
    }

    return V;
}
//=============================================================================

template<class isochrone_g>
VecDoub isochrone_grid<isochrone_g>::interp(double Z, double age, double M, std::vector<std::string> band){

    VecDoub V(2,-100.);

    int bot_Z,top_Z,bot_a,top_a;
    double dZ,mdZ,daD,mdaD,daU,mdaU,dm00,mdm00,dm10,mdm10,dm01,mdm01,dm11,mdm11;
    int bot_m00=0,top_m00=1,bot_m01=0,top_m01=1,bot_m10=0,top_m10=1,bot_m11=0,top_m11=1;

    if(Z<fehgrid.front()) Z=fehgrid.front();
    else if(Z>fehgrid.back()) Z=fehgrid.back();
    topbottom(fehgrid,Z,&bot_Z,&top_Z,"Isochrone Metal");
    dZ = (Z-fehgrid[bot_Z])/(fehgrid[top_Z]-fehgrid[bot_Z]); mdZ=1-dZ;

    if(age<agegrid.front()) age=agegrid.front();
    else if(age>agegrid.back()) age=agegrid.back();
    topbottom(agegrid,age,&bot_a,&top_a,"Isochrone Age");
    daU = (age-agegrid[bot_a])/(agegrid[top_a]-agegrid[bot_a]); mdaU=1-daU;
    daD = daU; mdaD = mdaU;

    top_m10=iso(top_Z,bot_a)->N()-1;
    if(M<iso(top_Z,bot_a)->initial_mass(bot_m10)) {dm10=0.;mdm10=1.;}
    else if(M>iso(top_Z,bot_a)->initial_mass(top_m10)) return V;
    else{
        topbottom(iso(top_Z,bot_a)->im(),M,&bot_m10,&top_m10,"Isochrone M");
        dm10 = (M-iso(top_Z,bot_a)->initial_mass(bot_m10))/(iso(top_Z,bot_a)->initial_mass(top_m10)-iso(top_Z,bot_a)->initial_mass(bot_m10));
        mdm10=1-dm10;
    }

    top_m11=iso(top_Z,top_a)->N()-1;
    if(M<iso(top_Z,top_a)->initial_mass(bot_m11)) {dm11=0.;mdm11=1.;}
    else if(M>iso(top_Z,top_a)->initial_mass(top_m11)){
        // do barycentric interp
        bot_m11=top_m11;
        VecDoub mm = {age,M};
        VecDoub mm1(2,0.);mm1[0]=agegrid[bot_a];mm1[1]=iso(top_Z,bot_a)->initial_mass(bot_m10);
        VecDoub mm2(2,0.);mm2[0]=agegrid[bot_a];mm1[1]=iso(top_Z,bot_a)->initial_mass(top_m10);
        VecDoub mm3(2,0.);mm3[0]=agegrid[top_a];mm1[1]=iso(top_Z,top_a)->initial_mass(bot_m11);
        VecDoub dm1 = mm1-mm2, dm2 = mm1-mm3;
        double afac = cross_product2D(dm1,dm2);
        mm1=mm1-mm;mm2 = mm2-mm;mm3=mm3-mm;
        mdm10 = cross_product2D(mm2,mm3)/afac;
        dm10 = cross_product2D(mm3,mm1)/afac;
        mdm11 = cross_product2D(mm1,mm2)/afac;
        dm11=0.;
        daU=1.;mdaU=1.;
    }
    else{
        topbottom(iso(top_Z,top_a)->im(),M,&bot_m11,&top_m11,"Isochrone M");
        dm11 = (M-iso(top_Z,top_a)->initial_mass(bot_m11))/(iso(top_Z,top_a)->initial_mass(top_m11)-iso(top_Z,top_a)->initial_mass(bot_m11));
        mdm11=1-dm11;
    }

    top_m00=iso(bot_Z,bot_a)->N()-1;
    if(M<iso(bot_Z,bot_a)->initial_mass(bot_m00)) {dm00=0.;mdm00=1.;}
    else if(M>iso(bot_Z,bot_a)->initial_mass(top_m00)) return V;
    else{
        topbottom(iso(bot_Z,bot_a)->im(),M,&bot_m00,&top_m00,"Isochrone M");
        dm00 = (M-iso(bot_Z,bot_a)->initial_mass(bot_m00))/(iso(bot_Z,bot_a)->initial_mass(top_m00)-iso(bot_Z,bot_a)->initial_mass(bot_m00));
        mdm00=1-dm00;
    }


    top_m01=iso(bot_Z,top_a)->N()-1;
    if(M<iso(bot_Z,top_a)->initial_mass(bot_m01)) {dm01=0.;mdm01=1.;}
    else if(M>iso(bot_Z,top_a)->initial_mass(top_m01)){
        // do barycentric interp
        bot_m01=top_m01;
        VecDoub mm = {age,M};
        VecDoub mm1(2,0.);mm1[0]=agegrid[bot_a];mm1[1]=iso(bot_Z,bot_a)->initial_mass(bot_m00);
        VecDoub mm2(2,0.);mm2[0]=agegrid[bot_a];mm1[1]=iso(bot_Z,bot_a)->initial_mass(top_m00);
        VecDoub mm3(2,0.);mm3[0]=agegrid[top_a];mm1[1]=iso(bot_Z,top_a)->initial_mass(bot_m01);
        VecDoub dm1 = mm1-mm2, dm2 = mm1-mm3;
        double afac = cross_product2D(dm1,dm2);
        mm1=mm1-mm;mm2 = mm2-mm;mm3=mm3-mm;
        mdm00 = cross_product2D(mm2,mm3)/afac;
        dm00 = cross_product2D(mm3,mm1)/afac;
        mdm01 = cross_product2D(mm1,mm2)/afac;
        dm01=0.;
        daD=1.;mdaD=1.;
    }
    else{
        topbottom(iso(bot_Z,top_a)->im(),M,&bot_m01,&top_m01,"Isochrone M");
        dm01 = (M-iso(bot_Z,top_a)->initial_mass(bot_m01))/(iso(bot_Z,top_a)->initial_mass(top_m01)-iso(bot_Z,top_a)->initial_mass(bot_m01));
        mdm01=1-dm01;
    }

    // Teff
    V[0] = dZ*(
               daU*(
                   dm11*iso(top_Z,top_a)->logTeff(top_m11)+
                  mdm11*iso(top_Z,top_a)->logTeff(bot_m11))
             +mdaU*(dm10*iso(top_Z,bot_a)->logTeff(top_m10)+
                  mdm10*iso(top_Z,bot_a)->logTeff(bot_m10)))
         +mdZ*(
               daD*(dm01*iso(bot_Z,top_a)->logTeff(top_m01)+
                  mdm01*iso(bot_Z,top_a)->logTeff(bot_m01))
             +mdaD*(dm00*iso(bot_Z,bot_a)->logTeff(top_m00)+
                  mdm00*iso(bot_Z,bot_a)->logTeff(bot_m00)));
    // logg
    V[1] = dZ*(
               daU*(
                   dm11*iso(top_Z,top_a)->logg(top_m11)+
                  mdm11*iso(top_Z,top_a)->logg(bot_m11))
             +mdaU*(dm10*iso(top_Z,bot_a)->logg(top_m10)+
                  mdm10*iso(top_Z,bot_a)->logg(bot_m10)))
         +mdZ*(
               daD*(
                   dm01*iso(bot_Z,top_a)->logg(top_m01)+
                  mdm01*iso(bot_Z,top_a)->logg(bot_m01))
             +mdaD*(dm00*iso(bot_Z,bot_a)->logg(top_m00)+
                  mdm00*iso(bot_Z,bot_a)->logg(bot_m00)));

    for(auto b: band)
        V.push_back(dZ*(daU*(dm11*iso(top_Z,top_a)->mag(top_m11,b)+
                      mdm11*iso(top_Z,top_a)->mag(bot_m11,b))
                 +mdaU*(dm10*iso(top_Z,bot_a)->mag(top_m10,b)+
                      mdm10*iso(top_Z,bot_a)->mag(bot_m10,b)))
             +mdZ*(daD*(dm01*iso(bot_Z,top_a)->mag(top_m01,b)+
                      mdm01*iso(bot_Z,top_a)->mag(bot_m01,b))
                 +mdaD*(dm00*iso(bot_Z,bot_a)->mag(top_m00,b)+
                      mdm00*iso(bot_Z,bot_a)->mag(bot_m00,b))));

    return V;
}
//=============================================================================

template<class isochrone_g>
void isochrone_grid<isochrone_g>::find_es_intercepts(double Z, double age, double es, int *bot_Z, int *top_Z, int *bot_a, int *top_a, int *bot_m00,int *top_m00,int *bot_m01,int *top_m01,int *bot_m10,int *top_m10,int *bot_m11,int *top_m11,double *dZ,double *mdZ,double *da,double *mda,double *dm00,double *mdm00,double *dm10,double *mdm10,double *dm01,double *mdm01,double *dm11,double *mdm11){

    if(Z<fehgrid.front()) Z=fehgrid.front();
    else if(Z>fehgrid.back()) Z=fehgrid.back();
    topbottom(fehgrid,Z,bot_Z,top_Z,"Isochrone Metal");
    *dZ = (Z-fehgrid[*bot_Z])/(fehgrid[*top_Z]-fehgrid[*bot_Z]); *mdZ=1-*dZ;

    if(age<agegrid.front()) age=agegrid.front();
    else if(age>agegrid.back()) age=agegrid.back();
    topbottom(agegrid,age,bot_a,top_a,"Isochrone Age");
    *da = (age-agegrid[*bot_a])/(agegrid[*top_a]-agegrid[*bot_a]); *mda=1-*da;


    if(es<iso(*top_Z,*top_a)->afraction()[0]){*bot_m11=0;*top_m11=1;}
    else topbottom(iso(*top_Z,*top_a)->afraction(),es,bot_m11,top_m11,"Isochrone es");
    *dm11 = (es-iso(*top_Z,*top_a)->age_fraction(*bot_m11))/(iso(*top_Z,*top_a)->age_fraction(*top_m11)-iso(*top_Z,*top_a)->age_fraction(*bot_m11));
    *mdm11=1.-*dm11;

    if(es<iso(*top_Z,*bot_a)->afraction()[0]){*bot_m10=0;*top_m10=1;}
    else topbottom(iso(*top_Z,*bot_a)->afraction(),es,bot_m10,top_m10,"Isochrone es");
    *dm10 = (es-iso(*top_Z,*bot_a)->age_fraction(*bot_m10))/(iso(*top_Z,*bot_a)->age_fraction(*top_m10)-iso(*top_Z,*bot_a)->age_fraction(*bot_m10));
    *mdm10=1.-*dm10;

    if(es<iso(*bot_Z,*top_a)->afraction()[0]){*bot_m01=0;*top_m01=1;}
    else topbottom(iso(*bot_Z,*top_a)->afraction(),es,bot_m01,top_m01,"Isochrone es");
    *dm01 = (es-iso(*bot_Z,*top_a)->age_fraction(*bot_m01))/(iso(*bot_Z,*top_a)->age_fraction(*top_m01)-iso(*bot_Z,*top_a)->age_fraction(*bot_m01));
    *mdm01=1.-*dm01;

    if(es<iso(*bot_Z,*bot_a)->afraction()[0]){*bot_m00=0;*top_m00=1;}
    else topbottom(iso(*bot_Z,*bot_a)->afraction(),es,bot_m00,top_m00,"Isochrone es");
    *dm00 = (es-iso(*bot_Z,*bot_a)->age_fraction(*bot_m00))/(iso(*bot_Z,*bot_a)->age_fraction(*top_m00)-iso(*bot_Z,*bot_a)->age_fraction(*bot_m00));
    *mdm00=1.-*dm00;
}
//=============================================================================

template<class isochrone_g>
VecDoub isochrone_grid<isochrone_g>::interp_es(double Z, double age, double M, std::vector<std::string> band){

    VecDoub V(2,-100.);

    int bot_Z,top_Z,bot_a,top_a;
    double dZ,mdZ,da,mda,dm00,mdm00,dm10,mdm10,dm01,mdm01,dm11,mdm11;
    int bot_m00=0,top_m00=1,bot_m01=0,top_m01=1,bot_m10=0,top_m10=1,bot_m11=0,top_m11=1;

    double es = age/max_age(Z,M); if(es>1.) return V;

    find_es_intercepts(Z,age,es,&bot_Z,&top_Z,&bot_a,&top_a,&bot_m00,&top_m00,&bot_m01,&top_m01,&bot_m10,&top_m10,&bot_m11,&top_m11,&dZ,&mdZ,&da,&mda,&dm00,&mdm00,&dm10,&mdm10,&dm01,&mdm01,&dm11,&mdm11);

    // Teff
    V[0] = dZ*(
               da*(
                   dm11*iso(top_Z,top_a)->logTeff(top_m11)+
                  mdm11*iso(top_Z,top_a)->logTeff(bot_m11))
             +mda*(dm10*iso(top_Z,bot_a)->logTeff(top_m10)+
                  mdm10*iso(top_Z,bot_a)->logTeff(bot_m10)))
         +mdZ*(
               da*(dm01*iso(bot_Z,top_a)->logTeff(top_m01)+
                  mdm01*iso(bot_Z,top_a)->logTeff(bot_m01))
             +mda*(dm00*iso(bot_Z,bot_a)->logTeff(top_m00)+
                  mdm00*iso(bot_Z,bot_a)->logTeff(bot_m00)));
    // logg
    V[1] = dZ*(
               da*(
                   dm11*iso(top_Z,top_a)->logg(top_m11)+
                  mdm11*iso(top_Z,top_a)->logg(bot_m11))
             +mda*(dm10*iso(top_Z,bot_a)->logg(top_m10)+
                  mdm10*iso(top_Z,bot_a)->logg(bot_m10)))
         +mdZ*(
               da*(
                   dm01*iso(bot_Z,top_a)->logg(top_m01)+
                  mdm01*iso(bot_Z,top_a)->logg(bot_m01))
             +mda*(dm00*iso(bot_Z,bot_a)->logg(top_m00)+
                  mdm00*iso(bot_Z,bot_a)->logg(bot_m00)));

    for(auto b: band)
        V.push_back(dZ*(da*(dm11*iso(top_Z,top_a)->mag(top_m11,b)+
                      mdm11*iso(top_Z,top_a)->mag(bot_m11,b))
                 +mda*(dm10*iso(top_Z,bot_a)->mag(top_m10,b)+
                      mdm10*iso(top_Z,bot_a)->mag(bot_m10,b)))
             +mdZ*(da*(dm01*iso(bot_Z,top_a)->mag(top_m01,b)+
                      mdm01*iso(bot_Z,top_a)->mag(bot_m01,b))
                 +mda*(dm00*iso(bot_Z,bot_a)->mag(top_m00,b)+
                      mdm00*iso(bot_Z,bot_a)->mag(bot_m00,b))));

    return V;
}
//=============================================================================

template<class isochrone_g>
double isochrone_grid<isochrone_g>::max_mass(double Z, double age){

    int bot_Z,top_Z,bot_a,top_a;
    double dZ,mdZ,da,mda;

    if(Z<fehgrid.front()) Z=fehgrid.front();
    else if(Z>fehgrid.back()) Z=fehgrid.back();
    topbottom(fehgrid,Z,&bot_Z,&top_Z,"MaxMass Isochrone Metal");
    dZ = (Z-fehgrid[bot_Z])/(fehgrid[top_Z]-fehgrid[bot_Z]); mdZ=1-dZ;

    if(age<agegrid.front()) age=agegrid.front();
    else if(age>agegrid.back()) age=agegrid.back();
    topbottom(agegrid,age,&bot_a,&top_a,"MaxMass Isochrone Age");
    da = (age-agegrid[bot_a])/(agegrid[top_a]-agegrid[bot_a]); mda=1-da;
    return dZ*(da*massmax[top_Z][top_a]+mda*massmax[top_Z][bot_a])
         +mdZ*(da*massmax[bot_Z][top_a]+mda*massmax[bot_Z][bot_a]);
}
//=============================================================================

template<class isochrone_g>
double isochrone_grid<isochrone_g>::max_age(double Z, double M){
    int bot_Z,top_Z,bot_M1,top_M1,bot_M2,top_M2;
    double dZ,mdZ,dM1,mdM1,dM2,mdM2;

    if(Z<fehgrid.front()) Z=fehgrid.front();
    else if(Z>fehgrid.back()) Z=fehgrid.back();
    topbottom(fehgrid,Z,&bot_Z,&top_Z,"MaxAge Isochrone Metal");
    dZ = (Z-fehgrid[bot_Z])/(fehgrid[top_Z]-fehgrid[bot_Z]); mdZ=1-dZ;

    if(M>massmax[bot_Z].front()){top_M1=1;bot_M1=0;}
    else if(M<massmax[bot_Z].back()){top_M1=agegrid.size()-1;bot_M1=top_M1-1;}
    else{
    topbottom(massmax[bot_Z],M,&bot_M1,&top_M1,"MaxAge Isochrone Mass");
    }
    dM1 = (log(M)-log(massmax[bot_Z][bot_M1]))/(log(massmax[bot_Z][top_M1])-log(massmax[bot_Z][bot_M1])); mdM1=1.-dM1;
    if(M>massmax[top_Z].front()){top_M2=1;bot_M2=0;}
    else if(M<massmax[top_Z].back()){top_M2=agegrid.size()-1;bot_M2=top_M2-1;}
    else{
    topbottom(massmax[top_Z],M,&bot_M2,&top_M2,"MaxAge Isochrone Mass");
    }
    dM2 = (log(M)-log(massmax[top_Z][bot_M2]))/(log(massmax[top_Z][top_M2])-log(massmax[top_Z][bot_M2])); mdM2=1.-dM2;
    return exp(mdZ*(dM1*log(agegrid[top_M1])+mdM1*log(agegrid[bot_M1]))
         +dZ*(dM2*log(agegrid[top_M2])+mdM2*log(agegrid[bot_M2])));
}
//=============================================================================

template<class isochrone_g>
double isochrone_grid<isochrone_g>::mass_from_es(double Z, double age, double es){
    int bot_Z,top_Z,bot_a,top_a;
    double dZ,mdZ,da,mda,dm00,mdm00,dm10,mdm10,dm01,mdm01,dm11,mdm11;
    int bot_m00=0,top_m00=1,bot_m01=0,top_m01=1,bot_m10=0,top_m10=1,bot_m11=0,top_m11=1;

    find_es_intercepts(Z,age,es,&bot_Z,&top_Z,&bot_a,&top_a,&bot_m00,&top_m00,&bot_m01,&top_m01,&bot_m10,&top_m10,&bot_m11,&top_m11,&dZ,&mdZ,&da,&mda,&dm00,&mdm00,&dm10,&mdm10,&dm01,&mdm01,&dm11,&mdm11);

    return dZ*(
               da*(
                   dm11*iso(top_Z,top_a)->initial_mass(top_m11)+
                  mdm11*iso(top_Z,top_a)->initial_mass(bot_m11))
             +mda*(dm10*iso(top_Z,bot_a)->initial_mass(top_m10)+
                  mdm10*iso(top_Z,bot_a)->initial_mass(bot_m10)))
         +mdZ*(
               da*(
                   dm01*iso(bot_Z,top_a)->initial_mass(top_m01)+
                  mdm01*iso(bot_Z,top_a)->initial_mass(bot_m01))
             +mda*(dm00*iso(bot_Z,bot_a)->initial_mass(top_m00)+
                  mdm00*iso(bot_Z,bot_a)->initial_mass(bot_m00)));
}
//=============================================================================

template<class isochrone_g>
double isochrone_grid<isochrone_g>::d_es_d_M(double Z, double age, double es){
    int bot_Z,top_Z,bot_a,top_a;
    double dZ,mdZ,da,mda,dm00,mdm00,dm10,mdm10,dm01,mdm01,dm11,mdm11;
    int bot_m00=0,top_m00=1,bot_m01=0,top_m01=1,bot_m10=0,top_m10=1,bot_m11=0,top_m11=1;

    find_es_intercepts(Z,age,es,&bot_Z,&top_Z,&bot_a,&top_a,&bot_m00,&top_m00,&bot_m01,&top_m01,&bot_m10,&top_m10,&bot_m11,&top_m11,&dZ,&mdZ,&da,&mda,&dm00,&mdm00,&dm10,&mdm10,&dm01,&mdm01,&dm11,&mdm11);

    double grad11 = (iso(top_Z,top_a)->age_fraction(top_m11)-iso(top_Z,top_a)->age_fraction(bot_m11))/(iso(top_Z,top_a)->initial_mass(top_m11)-iso(top_Z,top_a)->initial_mass(bot_m11));
    double grad10 = (iso(top_Z,bot_a)->age_fraction(top_m10)-iso(top_Z,bot_a)->age_fraction(bot_m10))/(iso(top_Z,bot_a)->initial_mass(top_m10)-iso(top_Z,bot_a)->initial_mass(bot_m10));
    double grad01 = (iso(bot_Z,top_a)->age_fraction(top_m01)-iso(bot_Z,top_a)->age_fraction(bot_m01))/(iso(bot_Z,top_a)->initial_mass(top_m01)-iso(bot_Z,top_a)->initial_mass(bot_m01));
    double grad00 = (iso(bot_Z,bot_a)->age_fraction(top_m00)-iso(bot_Z,bot_a)->age_fraction(bot_m00))/(iso(bot_Z,bot_a)->initial_mass(top_m00)-iso(bot_Z,bot_a)->initial_mass(bot_m00));

    return dZ*(da*grad11+mda*grad10)+mdZ*(da*grad01+mda*grad00);
}
//=============================================================================
template class isochrone_grid<isochrone_johnson>;
template class isochrone_grid<isochrone_padova>;
template class isochrone_grid<isochrone_dartmouth>;
//=============================================================================
