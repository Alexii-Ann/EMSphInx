#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "modality/laue/nml.hpp"

int main () {

emsphinx::laue::Namelist nml ;

//read default parameters from a single template namelist file
std::string nmlName = "RegularLaue.nml";                //specify the input/template namelist file
std::ifstream is(nmlName);                              //open the input file
std::string str((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());//read the input file into a string
std::string warning = nml.from_string(str);             //parse the string to a namelist

//list of bandwidths to produce namelist files for
std::vector<size_t> bandWidths = {54, 68, 88, 114, 124, 134, 140, 142, 144, 146, 152, 154, 158, 160, 180, 202, 204, 206, 234, 264, 284, 310, 338};

//loop over bandwidths producing files
for(size_t bw : bandWidths) {
    nml.bw = bw;
    std::ostringstream oss;//you need to #include <sstream> for this
    oss << "RegularLaue_" << bw << ".nml";
    std::ofstream os(oss.str());
    os << nml.to_string();
}


}

