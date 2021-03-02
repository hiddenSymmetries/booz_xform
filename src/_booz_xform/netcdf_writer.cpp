#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <netcdf.h>
#include "booz_xform.hpp"
#include "netcdf_writer.hpp"

// If we ever want to use a type other than double, we could change these next lines:
#define BOOZFLOAT NC_DOUBLE
#define nc_put_var_boozfloat nc_put_var_double

using namespace booz_xform;

booz_xform::NetCDFWriter::NetCDFWriter(std::string filename) {
  int retval;
  // If you want to use the netcdf4 format, "|NC_NETCDF4" must be included
  // in the file type, but such files cannot be read with
  // scipy.io.netcdf.
  if ((retval = nc_create(filename.c_str(), NC_CLOBBER, &ncid))) ERR(retval);
}

void booz_xform::NetCDFWriter::ERR(int e) {
  throw std::runtime_error(nc_strerror(e));
}

int booz_xform::NetCDFWriter::dim(std::string dimname, int val) {
  int dim_id, retval;
  if ((retval = nc_def_dim(ncid, dimname.c_str(), val, &dim_id)))
    ERR(retval);
  return dim_id;
}

/**
 * If an empty string is provided, no attribute is written.
 */
void booz_xform::NetCDFWriter::add_attribute(int var_id, std::string str, std::string units) {
  int retval;
  if (str.size() > 0) {
    if ((retval = nc_put_att_text(ncid, var_id, "description",
				  str.size(), str.c_str())))
      ERR(retval);
  }
  if (units.size() > 0) {
    if ((retval = nc_put_att_text(ncid, var_id, "units",
				  units.size(), units.c_str())))
      ERR(retval);
  }
}

void booz_xform::NetCDFWriter::put(std::string varname, int& val, std::string att, std::string units) {
  // Variant for scalar ints
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), NC_INT, 0, NULL, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(BOOZ_NC_INT);
  pointers.push_back((void*) &val);
  add_attribute(var_id, att, units);
}

void booz_xform::NetCDFWriter::put(std::string varname, boozfloat& val, std::string att, std::string units) {
  // Variant for scalar floats
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), BOOZFLOAT, 0, NULL, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(BOOZ_NC_FLOAT);
  pointers.push_back((void*) &val);
  add_attribute(var_id, att, units);
}

void booz_xform::NetCDFWriter::put(std::string varname, std::string& val, std::string att) {
  // Variant for strings. Note that unlike numerical quantities, strings have no units.
  int var_id, retval;

  // Make a name for the dimensions corresponding to the 
  int len = val.size();
  std::ostringstream converter;
  converter << "dim_" << std::setfill('0') << std::setw(5) << len;
  std::string dim_str = converter.str();

  // See if this dimension already exists:
  int dim_id;
  retval = nc_inq_dimid(ncid, dim_str.c_str(), &dim_id);
  if (retval != NC_NOERR) {
    // The dimension does not yet exist, so create it.
    if ((retval = nc_def_dim(ncid, dim_str.c_str(), len, &dim_id)))
      ERR(retval);
  }

  // Now that we have a dimension, define the string variable
  if ((retval = nc_def_var(ncid, varname.c_str(), NC_CHAR, 1, &dim_id, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(BOOZ_NC_STRING);
  pointers.push_back((void*) &val[0]);
  std::string units = "";
  add_attribute(var_id, att, units);
}

void booz_xform::NetCDFWriter::put(dim_id_type dim_id, std::string varname, IntVector& val, std::string att, std::string units) {
  // Variant for 1D int arrays
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), NC_INT, 1, &dim_id, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(BOOZ_NC_INT);
  pointers.push_back((void*) &val[0]);
  add_attribute(var_id, att, units);
}

void booz_xform::NetCDFWriter::put(dim_id_type dim_id, std::string varname, Vector& val, std::string att, std::string units) {
  // Variant for 1D float arrays
  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), BOOZFLOAT, 1, &dim_id, &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(BOOZ_NC_FLOAT);
  pointers.push_back((void*) &val[0]);
  add_attribute(var_id, att, units);
}

void booz_xform::NetCDFWriter::put(std::vector<dim_id_type> dim_id, std::string varname, boozfloat* pointer, std::string att, std::string units) {
  // Variant for ND float arrays for N > 1
  
  // NetCDF wants the order of the dimensions to be reversed compared to the booz_xform definitions.
  // Therefore we make a copy of the array of dimensions, and reverse the order of the copy.
  std::vector<dim_id_type> dim_id_reversed(dim_id);
  std::reverse(std::begin(dim_id_reversed), std::end(dim_id_reversed));

  int var_id, retval;
  if ((retval = nc_def_var(ncid, varname.c_str(), BOOZFLOAT, dim_id.size(), &dim_id_reversed[0], &var_id)))
    ERR(retval);
  var_ids.push_back(var_id);
  types.push_back(BOOZ_NC_FLOAT);
  pointers.push_back((void*) pointer);
  add_attribute(var_id, att, units);
}

void booz_xform::NetCDFWriter::write_and_close() {
  int retval;

  // End define mode. This tells netCDF we are done defining metadata.
  if ((retval = nc_enddef(ncid))) ERR(retval);

  // Write the data
  for (int j = 0; j < var_ids.size(); j++) {
    if (types[j] == BOOZ_NC_INT) {
      // ints
      if ((retval = nc_put_var_int(ncid, var_ids[j], (int*) pointers[j])))
	ERR(retval);
    } else if (types[j] == BOOZ_NC_STRING) {
      // strings
      if ((retval = nc_put_var_text(ncid, var_ids[j], (char*) pointers[j])))
	ERR(retval);
    } else {
      // floats/doubles:
      if ((retval = nc_put_var_boozfloat(ncid, var_ids[j], (boozfloat*) pointers[j])))
	ERR(retval);
    }
  }
  
  // Close the file
  if ((retval = nc_close(ncid))) ERR(retval);
}

