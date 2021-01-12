#include <iostream>
#include <stdexcept>
#include "booz_xform.hpp"
#include "netcdf_reader.hpp"

using namespace booz_xform;

/** Read in quantities for a Qsc object from a NetCDF file. This
    function is used for testing.
 */

void Booz_xform::read_boozmn(std::string filename) {
  if (verbose > 0) std::cout << "About to try reading netcdf file " << filename << std::endl;
  booz_xform::NetCDFReader nc(filename);

  nc.get("mboz_b", mboz);
  nc.get("nboz_b", nboz);

  if (verbose) std::cout << "Read mboz=" << mboz << " , nboz=" << nboz << std::endl;
  
  nc.close();
  
}
