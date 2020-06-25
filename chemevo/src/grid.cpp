#include "grid.h"
//=============================================================================
// Constructors
//=============================================================================
Grid::Grid(unsigned NR, unsigned NT,double Rmin, double Rmax, double MaxAge){
	radial_grid = create_range(Rmin,Rmax,NR);
	time_grid = create_range(0.,MaxAge,NT);
	age_grid = create_range(MaxAge,0.,NT);
	grid = std::vector<VecDoub>(NR,VecDoub(NT,0.));
}
Grid::Grid(ModelParameters params){
	MaxAge = extract_param(params.parameters["fundamentals"],"GalaxyAge", 13.7);
	double Rmin = extract_param(params.parameters["grids"],
	                            "MinimumRadius", 0.1);
	double Rmax = extract_param(params.parameters["grids"],
	                            "MaximumRadius", 15.);
	int NR = extract_param(params.parameters["grids"], "RadialGridPoints", 50);
	int NT = extract_param(params.parameters["grids"], "AgeGridPoints", 50);
	radial_grid = create_range(Rmin,Rmax,NR);
	if(NR==1)
		radial_grid=VecDoub(1,.5*(Rmax+Rmin));
	bool LAGrid = extract_param(params.parameters["grids"],
	                            "LogAgeGrid", true);
	if(LAGrid){
		time_grid = create_log_range(3e-3,MaxAge,NT);
		time_grid[0]=0.;
		}
	else
		time_grid = create_range(0.,MaxAge,NT);
	age_grid = time_grid;
	for(int i=0; i<NT; ++i) age_grid[i]=MaxAge-age_grid[i];
	grid = std::vector<VecDoub>(NR,VecDoub(NT,0.));
}
//=============================================================================
// Getter
//=============================================================================
double Grid::operator()(double r, double t, bool loginterp, bool extrapolate){
	if(t<time_grid.front())t=time_grid.front();
	if(t>time_grid.back())t=time_grid.back();
	if(r<radial_grid.front() and radial_grid.size()>1){
		if(extrapolate) return log_extrapolate_low(r,t);
		else r=radial_grid.front();
	}
	if(r>radial_grid.back() and radial_grid.size()>1){
		if(extrapolate) return log_extrapolate_high(r,t);
		else r=radial_grid.back();
	}
	int bott, topt, botr, topr;
	if(radial_grid.size()>1) topbottom(radial_grid,r,&botr,&topr);
        else {topr=0;botr=0;}
	topbottom(time_grid,t,&bott,&topt);
	double y1 = grid[botr][bott]+(t-time_grid[bott])/(time_grid[topt]-time_grid[bott])*(grid[botr][topt]-grid[botr][bott]);
	double y2 = grid[topr][bott]+(t-time_grid[bott])/(time_grid[topt]-time_grid[bott])*(grid[topr][topt]-grid[topr][bott]);
        if(loginterp){
		y1=log(y1);
		y2=log(y2);
	}
	double rslt=y1;
        if(radial_grid.size()>1)
            rslt += (r-radial_grid[botr])/(radial_grid[topr]-radial_grid[botr])*(y2-y1);
	if(loginterp)
		return exp(rslt);
        return rslt;
}
// Extrapolation
double Grid::log_extrapolate_low(double R, double t){
	double Rlow = radial_grid[0], Rhigh = radial_grid[1];
	return exp(log((*this)(Rlow,t))+(log((*this)(Rhigh,t))-log((*this)(Rlow,t)))/(Rhigh-Rlow)*(R-Rlow));
}
double Grid::log_extrapolate_low(double R, unsigned t){
	double Rlow = radial_grid[0], Rhigh = radial_grid[1];
	return exp(log(grid[0][t])+(log(grid[1][t])-log(grid[0][t]))/(Rhigh-Rlow)*(R-Rlow));
}
double Grid::log_extrapolate_high(double R, double t){
	unsigned NR=radial_grid.size();
	double Rlow = radial_grid[NR-2], Rhigh = radial_grid[NR-1];
	return exp(log((*this)(Rhigh,t))+(log((*this)(Rhigh,t))-log((*this)(Rlow,t)))/(Rhigh-Rlow)*(R-Rhigh));
}
double Grid::log_extrapolate_high(double R, unsigned t){
	unsigned NR=radial_grid.size();
	double Rlow = radial_grid[NR-2], Rhigh = radial_grid[NR-1];
	return exp(log(grid[NR-1][t])+(log(grid[NR-1][t])-log(grid[NR-2][t]))/(Rhigh-Rlow)*(R-Rhigh));
}

double Grid::t_gradient(double r, double t){
	if(r<radial_grid.front())r=radial_grid.front();
	if(r>radial_grid.back())r=radial_grid.back();
	if(t<time_grid.front())t=time_grid.front();
	if(t>time_grid.back())t=time_grid.back();
	int bott, topt, botr, topr;
	topbottom(radial_grid,r,&botr,&topr);
	topbottom(time_grid,t,&bott,&topt);
	double y1 = (grid[botr][topt]-grid[botr][bott])/(time_grid[topt]-time_grid[bott]);
	double y2 = (grid[topr][topt]-grid[topr][bott])/(time_grid[topt]-time_grid[bott]);
	return y1+(r-radial_grid[botr])/(radial_grid[topr]-radial_grid[botr])*(y2-y1);
}
double Grid::R_gradient(double r, double t){
        if(radial_grid.size()==1) return 0.;
	if(r<radial_grid.front())r=radial_grid.front();
	if(r>radial_grid.back())r=radial_grid.back();
	if(t<time_grid.front())t=time_grid.front();
	if(t>time_grid.back())t=time_grid.back();
	int bott, topt, botr, topr;
	topbottom(radial_grid,r,&botr,&topr);
	topbottom(time_grid,t,&bott,&topt);
	double y1 = (grid[topr][bott]-grid[botr][bott])/(radial_grid[topr]-radial_grid[botr]);
	double y2 = (grid[topr][topt]-grid[botr][topt])/(radial_grid[topr]-radial_grid[botr]);
	return y1+(t-time_grid[bott])/(time_grid[topt]-time_grid[bott])*(y2-y1);
}
//=============================================================================
// Setters
//=============================================================================
void Grid::set_fixed_r_t_const(double c){
	for(unsigned i=0;i<radial_grid.size();++i)
	for(unsigned j=0;j<time_grid.size();++j)
		grid[i][j]=c;
}
void Grid::set_fixed_t_const(double c,unsigned t){
	assert(t<time_grid.size());
	for(unsigned i=0;i<radial_grid.size();++i)
		grid[i][t]=c;
}
void Grid::set_fixed_t(VecDoub rad,unsigned t){
	assert(rad.size()==radial_grid.size());
	assert(t<time_grid.size());
	for(unsigned i=0;i<rad.size();++i)
		grid[i][t]=rad[i];
}
void Grid::set_fixed_r(VecDoub tim,unsigned r){
	assert(tim.size()==time_grid.size());
	assert(r<radial_grid.size());
	for(unsigned i=0;i<tim.size();++i)
		grid[r][i]=tim[i];
}
void Grid::set(double x, unsigned r, unsigned t){
	assert(t<time_grid.size());
	assert(r<radial_grid.size());
	grid[r][t]=x;
}
void Grid::add_time(unsigned nt, double t){
	for(unsigned nr=0;nr<radial_grid.size();++nr)
		grid[nr].insert(grid[nr].begin()+nt,0.);
	time_grid.insert(time_grid.begin()+nt,t);
	age_grid.insert(age_grid.begin()+nt,MaxAge-t);
}
//=============================================================================
// Grid spacing
double Grid::Rdown(unsigned nR){
	double Rd=(nR>0?radial_grid[nR-1]:2.*radial_grid[0]-radial_grid[1]);
	if(Rd<0.) Rd=radial_grid[0]*0.5;
	return Rd;
}
double Grid::Rup(unsigned nR){
	auto NR = radial_grid.size();
	double Ru=(nR<NR-1?radial_grid[nR+1]:2.*radial_grid[nR]-radial_grid[nR-1]);
	return Ru;
}
double Grid::annulus_area(int nR){
	double R, Rd, Ru;
	if(nR==radial_grid.size()){
		R = Rup(nR-1);
		Rd=radial_grid[nR-1];
		Ru = 2.*R-Rd;
	}
	else if(nR<0){
		Ru = radial_grid[0];
		R = Rdown(0);
		Rd = 2.*R-Ru;
		if(Rd<0.) Rd=0.;
	}
	else{
		Ru=Rup(nR);
		Rd=Rdown(nR);
		R =radial_grid[nR];
	}
	return .25*PI*(Ru*Ru-Rd*Rd+2.*R*(Ru-Rd));
}
double Grid::annulus_width(int nR){
	double R, Rd, Ru;
	if(nR==radial_grid.size()){
		R = Rup(nR-1);
		Rd=radial_grid[nR-1];
		Ru = 2.*R-Rd;
	}
	else if(nR<0){
		Ru = radial_grid[0];
		R = Rdown(0);
		Rd = 2.*R-Ru;
		if(Rd<0.) Rd=0.;
	}
	else{
		Ru=Rup(nR);
		Rd=Rdown(nR);
		R =radial_grid[nR];
	}
	return .5*(Ru-Rd);
}
//=============================================================================
// Printing
//=============================================================================
void Grid::print_grid(std::ostream &out){
	auto nR=0;
	out<<"grid ";
	for(auto t:time_grid)
		out<<t<<" ";
	out<<std::endl;
	for(auto g:grid){
		out<<radial_grid[nR]<<" ";
		for(auto t: g)
			out<<t<<" ";
		out<<std::endl;
		nR++;
	}
}

void Grid::print_grid_fixed_radius(std::ostream &out,double R){
	for(auto t:time_grid)
		out<<(*this)(R,t)<<" ";
	out<<std::endl;
}
void Grid::print_grid_fixed_time(std::ostream &out,double t){
	for(auto r:radial_grid)
		out<<(*this)(r,t)<<" ";
	out<<std::endl;
}
void Grid::print_grid_now(std::ostream &out){
	for(auto nR=0u;nR<radial_grid.size();++nR)
		out<<grid[nR].back()<<" ";
	out<<std::endl;
}
//=============================================================================
// Writing to HDF5
//=============================================================================
using namespace H5;

void Grid::write_radial_grid_hdf5(H5File &fout){
	hdf5_write_1D_vector(fout,radial_grid,"R");
}
void Grid::write_time_grid_hdf5(H5File &fout){
	hdf5_write_1D_vector(fout,time_grid,"t");
}
void Grid::write_hdf5(H5File &fout,std::string name){
	hdf5_write_2D(fout,grid,name);
}
//=============================================================================
