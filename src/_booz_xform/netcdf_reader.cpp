#include <iostream>
#include <stdexcept>
#include <netcdf.h>
#include "vector_matrix.hpp"
#include "netcdf_reader.hpp"

using namespace booz_xform;

booz_xform::NetCDFReader::NetCDFReader(std::string filename) {
  int retval;
  if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
    ERR(retval, filename);
  // Get the number of dimensions, variables, etc:
  if ((retval = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
    ERR(retval, filename);
}

void booz_xform::NetCDFReader::ERR(int e, std::string details) {
  std::cout << "NetCDF Error! " << nc_strerror(e)
            << " [" << details << "]" << std::endl;
  throw std::runtime_error(nc_strerror(e));
}

size_t booz_xform::NetCDFReader::getdim(std::string dimname) {
  int dim_id, retval;
  size_t dimval;
  if ((retval = nc_inq_dimid(ncid, dimname.c_str(), &dim_id)))
    ERR(retval, dimname);
  if ((retval = nc_inq_dimlen(ncid, dim_id, &dimval)))
    ERR(retval, dimname);
  return dimval;
}

void booz_xform::NetCDFReader::get(std::string varname, int& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
    ERR(retval, varname);
  if ((retval = nc_get_var_int(ncid, var_id, &var)))
    ERR(retval, varname);
}

void booz_xform::NetCDFReader::get(std::string varname, boozfloat& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
    ERR(retval, varname);
  if ((retval = nc_get_var_double(ncid, var_id, &var)))
    ERR(retval, varname);
}

/** This subroutine presently assumes the Vector has already been
    sized to the correct dimension. If not, there may be a seg fault.
 */
void booz_xform::NetCDFReader::get(std::string varname, Vector& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
    ERR(retval, varname);
  if ((retval = nc_get_var_double(ncid, var_id, &var[0])))
    ERR(retval, varname);
}

/** This subroutine presently assumes the IntVector has already been
    sized to the correct dimension. If not, there may be a seg fault.
 */
void booz_xform::NetCDFReader::get(std::string varname, IntVector& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
    ERR(retval, varname);
  if ((retval = nc_get_var_int(ncid, var_id, &var[0])))
    ERR(retval, varname);
}

/** This subroutine presently assumes the matrix has already been
    sized to the correct dimension. If not, there may be a seg fault.
 */
void booz_xform::NetCDFReader::get(std::string varname, Matrix& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
    ERR(retval, varname);
  if ((retval = nc_get_var_double(ncid, var_id, &var(0,0))))
    ERR(retval, varname);
}

void booz_xform::NetCDFReader::close() {
  int retval;
  if ((retval = nc_close(ncid))) ERR(retval, "close()");
}
