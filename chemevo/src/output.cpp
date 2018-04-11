#include "model.h"

void Model::print_abundance_grid(std::ostream &out,Element E){
	mass_fraction[elements_r[E]].print_grid(out);
}
void Model::print_abundance_grid_fixed_radius(std::ostream &out,Element E,double R){
	mass_fraction[elements_r[E]].print_grid_fixed_radius(out,R);
}
void Model::print_abundance_grid_fixed_time(std::ostream &out,Element E,double t){
	mass_fraction[elements_r[E]].print_grid_fixed_time(out,t);
}
void Model::print_abundance_grid_now(std::ostream &out,Element E){
	mass_fraction[elements_r[E]].print_grid_now(out);
}

void Model::print_gasmass_grid(std::ostream &out){
	gas_mass->print_grid(out);
}
void Model::print_gasmass_grid_fixed_radius(std::ostream &out,double R){
	gas_mass->print_grid_fixed_radius(out,R);
}
void Model::print_gasmass_grid_fixed_time(std::ostream &out,double t){
	gas_mass->print_grid_fixed_time(out,t);
}
void Model::print_gasmass_grid_now(std::ostream &out){
	gas_mass->print_grid_now(out);
}

void Model::print_metallicity_grid(std::ostream &out){
	metallicity->print_grid(out);
}
void Model::print_metallicity_grid_fixed_radius(std::ostream &out,double R){
	metallicity->print_grid_fixed_radius(out,R);
}
void Model::print_metallicity_grid_fixed_time(std::ostream &out,double t){
	metallicity->print_grid_fixed_time(out,t);
}
void Model::print_metallicity_grid_now(std::ostream &out){
	metallicity->print_grid_now(out);
}
