#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include "../data/toml.hpp"

using namespace std;

DataFile::DataFile(std::string file_name)
: _file_name(file_name)
{
   // Lecture du fichier de données
   auto config = toml::parse(file_name);

   // Other
   const auto& parameter = toml::find(config, "parameter");
   this->_cas = toml::find<int>(parameter, "cas");
   this->_xmin = toml::find<double>(parameter, "xmin");
   this->_xmax = toml::find<double>(parameter, "xmax");
   this->_ymin = toml::find<double>(parameter, "ymin");
   this->_ymax = toml::find<double>(parameter, "ymax");
   this->_D = toml::find<double>(parameter, "D");
   this->_Tf = toml::find<double>(parameter, "Tf");
   this->_Nx = toml::find<int>(parameter, "Nx");
   this->_Ny = toml::find<int>(parameter, "Ny");
   this->_dt = toml::find<double>(parameter, "dt");
   this->_solver = toml::find<std::string>(parameter, "Solver");
   this->_r = toml::find<int>(parameter, "r");
   this->_alpha = toml::find<double>(parameter, "alpha");
   this->_beta = toml::find<double>(parameter, "beta");
}

#define _DATA_FILE_CPP
#endif
