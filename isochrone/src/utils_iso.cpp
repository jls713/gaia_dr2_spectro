#include "utils_iso.h"
//=============================================================================
void GetFilesInDirectory(std::vector<std::string> &out, const std::string &directory){
	/* Returns a list of files in a directory (except the ones that begin with a dot) */
    DIR *dir;
    struct dirent *ent;
    struct stat st;

    dir = opendir(directory.c_str());
    while ((ent = readdir(dir)) != NULL) {
        const std::string file_name = ent->d_name;
        const std::string full_file_name = directory + "/" + file_name;
        if (file_name[0] == '.') continue;
        if (stat(full_file_name.c_str(), &st) == -1) continue;
        const bool is_directory = (st.st_mode & S_IFDIR) != 0;
        if (is_directory) continue;
        out.push_back(full_file_name);
    }
    closedir(dir);
}
//=============================================================================
double KroupaIMF_default(double M){
    // This is Kroupa, Tout & Gilmore (1993) -- differs in high mass end from 
    // Kroupa (2001) by -2.7 -> -2.3
	if(M>=0.08 and M<0.5) return 0.035*pow(M,-1.3);
	else if(M>=0.5 and M<1.0) return 0.019*pow(M,-2.2);
	else  if(M>=1.0)  return 0.019*pow(M,-2.7);
    else return 0.;
}
//=============================================================================
double I2MASS(double J, double K){
    return 2*J-K+0.2*exp(((J-K)-1.2)/0.2)+0.12;
}
//=============================================================================
double cross_product2D(const VecDoub &a, const VecDoub &b){
    return a[0]*b[1]-b[0]*a[1];
}
//=============================================================================
double VpMag(double B, double V){
    return V+0.57*(B-V);
}
//=============================================================================
// From CASU pages
VecDoub JKH_vista_from_2MASS(double J2, double K2, double H2){
    return
    {
     J2-0.031*(J2-K2),
     K2-0.006*(J2-K2),
     H2+0.015*(J2-K2)
    };
}
//=============================================================================
// From Koen (2007)
VecDoub JKH_2MASS(double J, double K, double H){
    double K2MASS = K-0.050*(J-H)*(J-H);
    return
    {
     K2MASS + 0.859*(J-K) + 0.156*(J-H)*(J-H),
     K2MASS,
     K2MASS + 0.046 + 0.944*(H-K)
    };
}
//=============================================================================
json parameters(void){
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
std::vector<std::string> Split(const std::string& strIn,
				const std::string& token) {
   std::string str = strIn;
   std::vector<std::string>result;
    while(str.size()){
        int index = str.find(token);
        if(index!=std::string::npos){
            result.push_back(str.substr(0,index));
            str = str.substr(index+token.size());
            if(str.size()==0)result.push_back(str);
        }else{
            result.push_back(str);
            str = "";
        }
    }
    return result;
}
//=============================================================================
