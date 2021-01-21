#ifndef BOOZ_NETCDF_WRITER_H
#define BOOZ_NETCDF_WRITER_H

#include <vector>
#include <netcdf.h>
#include "booz_xform.hpp"

/** Here we use the C interface to NetCDF, not the C++ interface. One
    reason for this is to simplify the build system and avoid a
    potential problem for users, since the C++ interface does not
    automatically come with NetCDF. Second, I never managed to get the
    build system to work with the C++ interface on my laptop. Third,
    NetCDF's C++ interface is rather clunky, requiring multiple lines
    of code to add each variable, so I wanted to write a wrapper class
    to simplify this anyway.
 */

// using namespace booz_xform;

namespace booz_xform {
  
  typedef int dim_id_type;

  /** A class to streamline the process of writing a NetCDF file.
   */
  class NetCDFWriter {
  private:
    int ncid;
    std::vector<int> var_ids;
    std::vector<void*> pointers;
    enum {BOOZ_NC_INT, BOOZ_NC_FLOAT, BOOZ_NC_STRING};
    std::vector<int> types;
    static void ERR(int);
    
  public:
    NetCDFWriter(std::string);
    dim_id_type dim(std::string, int);
    void add_attribute(int, std::string, std::string);
    
    // Scalars:
    void put(std::string, int&, std::string, std::string);
    void put(std::string, boozfloat&, std::string, std::string);
    void put(std::string, std::string&, std::string);
    // 1D vectors:
    void put(dim_id_type, std::string, IntVector&, std::string, std::string);
    void put(dim_id_type, std::string, Vector&, std::string, std::string);
    // ND vectors for N > 1:
    void put(std::vector<dim_id_type>, std::string, boozfloat*, std::string, std::string);
    
    void write_and_close();
  };
}

#endif
