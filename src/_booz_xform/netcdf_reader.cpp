#include <iostream>
#include <stdexcept>
#include <netcdf.h>
#include "vector_matrix.hpp"
#include "netcdf_reader.hpp"

using namespace booz_xform;

booz_xform::NetCDFReader::NetCDFReader(std::string filename) {
  int retval;
  if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
    ERR(retval);
  // Get the number of dimensions, variables, etc:
  if ((retval = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
      ERR(retval);
}

void booz_xform::NetCDFReader::ERR(int e) {
  std::cout << "NetCDF Error! " << nc_strerror(e) << std::endl;
  throw std::runtime_error(nc_strerror(e));
}

void booz_xform::NetCDFReader::get(std::string varname, int& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_int(ncid, var_id, &var)))
      ERR(retval);
}

void booz_xform::NetCDFReader::get(std::string varname, boozfloat& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_double(ncid, var_id, &var)))
      ERR(retval);
}

/** This subroutine presently assumes the Vector has already been
    sized to the correct dimension. If not, there may be a seg fault.
 */
void booz_xform::NetCDFReader::get(std::string varname, Vector& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_double(ncid, var_id, &var[0])))
      ERR(retval);
}

/** This subroutine presently assumes the IntVector has already been
    sized to the correct dimension. If not, there may be a seg fault.
 */
void booz_xform::NetCDFReader::get(std::string varname, IntVector& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_int(ncid, var_id, &var[0])))
      ERR(retval);
}

/** This subroutine presently assumes the matrix has already been
    sized to the correct dimension. If not, there may be a seg fault.
 */
void booz_xform::NetCDFReader::get(std::string varname, Matrix& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_double(ncid, var_id, &var(0,0))))
      ERR(retval);
}

void booz_xform::NetCDFReader::close() {
  int retval;
  if ((retval = nc_close(ncid))) ERR(retval);
}
