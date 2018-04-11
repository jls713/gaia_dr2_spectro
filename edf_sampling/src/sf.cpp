#include "sf.h"

Binned_selection_function::Binned_selection_function(unsigned N_mag, double magmin, double magmax, std::string band, int NSIDE, std::string filename):N_mag(N_mag),NSIDE(NSIDE),band(band),filename(filename){
    Healpix_Ordering_Scheme scheme = RING;//NEST;
    for(unsigned i=0;i<N_mag;i++){
        HP_Map_grid.push_back(Healpix_Map<double>());
        HP_Map_grid[i].SetNside(NSIDE,scheme);
        read_Healpix_map_from_fits<double>(filename,HP_Map_grid[i],i+1);
    }
    npix=12*NSIDE*NSIDE;
    // These are bin edges
    mag_edges = create_range<double>(magmin,magmax,N_mag+1);
    for(unsigned i=0;i<N_mag;++i)
      mag.push_back((mag_edges[i]+mag_edges[i+1])*.5);
    std::cerr<<"Selection function loaded"<<std::endl;
}

double Binned_selection_function::evaluate(double l, double b, VecDoub m,bool interp){
  double M = m[0];
  if(b>PI/2. or b<-PI/2.) return 0.;
  b=PI/2.-b;
  int bot_M, top_M;
  // bot_I and top_I are bin edges. Bin value stored in bot_I element
  if(M<mag_edges.front() or M>mag_edges.back()) return 0.;
  if(M<mag.front()){
    if(interp) return HP_Map_grid[0].interpolated_value(pointing(b,l));
    else return HP_Map_grid[0][HP_Map_grid[0].ang2pix(pointing(b,l))];
  }
  if(M>mag.back()){
    if(interp) return HP_Map_grid[N_mag-1].interpolated_value(pointing(b,l));
    else return HP_Map_grid[N_mag-1][HP_Map_grid[N_mag-1].ang2pix(pointing(b,l))];
  }
  topbottom<double>(mag,M,&bot_M,&top_M,"SF M");
  double yu,yd;
  if(interp){
    yu = HP_Map_grid[top_M].interpolated_value(pointing(b,l));
    yd = HP_Map_grid[bot_M].interpolated_value(pointing(b,l));
  }
  else{
    yu = HP_Map_grid[top_M][HP_Map_grid[top_M].ang2pix(pointing(b,l))];
    yd = HP_Map_grid[bot_M][HP_Map_grid[bot_M].ang2pix(pointing(b,l))];
  }
  return yd+(M-mag[bot_M])/(mag[top_M]-mag[bot_M])*(yu-yd);
  // return HP_Map_grid[bot_M][HP_Map_grid[bot_M].ang2pix(pointing(b,l))];
  // return HP_Map_grid[bot_M].interpolated_value(pointing(b,l));
}

double Binned_selection_function::evaluate_at_pixel(double M, int pix, double *l, double *b){
  int bot_M, top_M;
  // bot_M and top_M are bin edges. Bin value stored in bot_M element
  if(M<mag_edges.front() or M>mag_edges.back()) return 0.;
  topbottom<double>(mag_edges,M,&bot_M,&top_M,"SF M");
  pointing p = HP_Map_grid[bot_M].pix2ang(pix);
  *b=PI/2.-p.theta;*l=p.phi;
  return HP_Map_grid[bot_M][pix];
}

double Binned_selection_function::pixel_to_lb(int pix, double *l, double *b){
  pointing p = HP_Map_grid[0].pix2ang(pix);
  *b=PI/2.-p.theta;*l=p.phi;
}

// int main(){
//     Binned_selection_function RSF;
//     std::cout<<RSF.evaluate(-0.872,-1.5,11.7);
// }
