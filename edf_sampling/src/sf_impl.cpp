#include "sf_impl.h"

std::string jason_dir(void){
    std::string s=parameters()["selection_function_folder"];
    return s+"TGAS/";
}
std::string gus_dir(void){
    std::string s=parameters()["selection_function_folder"];
    return s+"gus/";}
std::string rave_dir(void){
    std::string s=parameters()["selection_function_folder"];
    return s+"RAVE/";}

std::map<std::string, Selection_function*(*)(void)> sf_types={
  {"RAVE",&createSFInstance<RAVE_selection_function>},
  {"TGAS",&createSFInstance<TGAS_selection_function>},
  {"TGAS_HQ",&createSFInstance<TGAS_HQ_selection_function>},
  {"Tycho2",&createSFInstance<Tycho2_selection_function>},
  {"Tycho2_8",&createSFInstance<Tycho2_selection_function8>},
  {"Tycho2_32",&createSFInstance<Tycho2_selection_function32>},
  {"Tycho2_RAVE",&createSFInstance<Tycho2RAVE_selection_function>},
  {"Tycho2_RAVE_8",&createSFInstance<Tycho2RAVE_selection_function8>},
  {"Tycho2_RAVE_32",&createSFInstance<Tycho2RAVE_selection_function32>},
  {"Tycho2_RAVE_I",&createSFInstance<Tycho2RAVE_selection_function_I>},
  {"Tycho2_RAVE_8_I",&createSFInstance<Tycho2RAVE_selection_function8_I>},
  {"Tycho2_RAVE_32_I",&createSFInstance<Tycho2RAVE_selection_function32_I>}};
