#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "booz_xform.hpp"

//using namespace booz_xform;

int booz_xform::driver(int argc, char* argv[]) {
  int j;

  if (argc != 2) {
    std::cout << "Usage:  xbooz_xform <inputfile>" << std::endl;
    std::cout << std::endl;
    std::cout << "Where <inputfile> is a text file with 2 or 3 lines, formatted as follows:" << std::endl;
    std::cout << std::endl;
    std::cout << "<mboz> <nboz>" << std::endl;
    std::cout << "<extension>" << std::endl;
    std::cout << "<surfaces>" << std::endl;
    std::cout << std::endl;
    std::cout << "Here, " << std::endl;
    std::cout << std::endl;
    std::cout << "<mboz> is the maximum poloidal Fourier mode number to include in the Boozer representation." << std::endl;
    std::cout << std::endl;
    std::cout << "<nboz> is the maximum toroidal Fourier mode number to include in the Boozer representation." << std::endl;
    std::cout << std::endl;
    std::cout << "<extension> is the VMEC extension from the VMEC input file to process." << std::endl;
    std::cout << "For instance, to process the file wout_li383_1.4m.nc, <extension> should be li383_1.4m" << std::endl;
    std::cout << "Note that the wout file to load must be in the current working directory." << std::endl;
    std::cout << std::endl;
    std::cout << "<surfaces> is a list of integers giving the surfaces of the VMEC input file to" << std::endl;
    std::cout << "transform. Note that the transformation is only performed on half-grid surfaces." << std::endl;
    std::cout << "The first half-grid surface has index 2. The outermost available surface has index NS," << std::endl;
    std::cout << "where NS is the VMEC input parameter. You can omit <surfaces> from the input file," << std::endl;
    std::cout << "in which case the transformation will be performed on all half-grid surfaces." << std::endl;
    
    return 0;
  }
  
  std::cout << "This is xbooz_xform." << std::endl;

  // Read the input file.
  
  std::ifstream file;
  file.open(argv[1]);
  if (!file.is_open()) throw std::runtime_error(std::string("Unable to open input file ") + argv[1]);

  int mboz_in, nboz_in;
  file >> mboz_in >> nboz_in;
  std::cout << "Read mboz = " << mboz_in << ", nboz = " << nboz_in << std::endl;
  if (file.fail()) throw std::runtime_error("Unable to read mboz and nboz");

  std::string extension;
  file >> extension;
  std::cout << "Read extension = " << extension << std::endl;
  if (file.fail()) throw std::runtime_error("Unable to read extension");

  std::vector<int> jlist;
  int val;
  while (true) {
    file >> val;
    if (file.fail()) {
      std::cout << "Unable to read more." << std::endl;
      break;
    }
    std::cout << "Read value " << val << std::endl;
    jlist.push_back(val);
  }
  std::cout << "Read jlist =";
  for (j = 0; j < jlist.size(); j++) std::cout << " " << jlist[j];
  std::cout << std::endl;  
  
  file.close();

  // Done reading input file. Now set up the calculation.
  
  booz_xform::Booz_xform booz;
  booz.read_wout("wout_" + extension + ".nc");

  booz.mboz = std::max(booz.mboz, mboz_in);
  booz.nboz = std::max(booz.nboz, nboz_in);

  // If no surfaces are specified, we do not need to modify the
  // default jlist initialized by read_netcdf().
  if (jlist.size() > 0) {
    std::sort(jlist.begin(), jlist.end());
    booz.jlist.resize(jlist.size());
    for (j = 0; j < jlist.size(); j++) booz.jlist[j] = jlist[j];
  } else {
    std::cout << "No jlist specified, so including all half-grid surfaces." << std::endl;
  }
  std::cout << "About to run transformation with jlist =";
  for (j = 0; j < booz.jlist.size(); j++) std::cout << " " << booz.jlist[j];
  std::cout << std::endl;
  
  // Run the main calculation:
  booz.run();

  // Save results:
  booz.write_boozmn("boozmn_" + extension + ".nc");
  
  std::cout << "Good bye." << std::endl;

  return 0;
}
