#ifndef BOOZ_NETCDF_READER_H
#define BOOZ_NETCDF_READER_H

#include "vector_matrix.hpp"

namespace booz_xform {
  
  /** A class to streamline the process of reading a NetCDF file.
   */
  class NetCDFReader {
  private:
    int ncid, ndims, nvars, ngatts, unlimdimid;
    static void ERR(int, std::string);
    
  public:
    NetCDFReader(std::string);

    size_t getdim(std::string);
    
    // Scalars:
    void get(std::string, int&);
    void get(std::string, boozfloat&);
    // Vectors
    void get(std::string, Vector&);
    void get(std::string, IntVector&);
    // Matrices
    void get(std::string, Matrix&);
	     
    void close();
  };
}

#endif
