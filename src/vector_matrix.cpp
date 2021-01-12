#include <iostream>
#include <valarray>
#include <cassert>
#include "vector_matrix.hpp"

// Representation of BLAS and LAPACK routines we need:
extern "C" {
  // Matrix-vector multiply:
  void dgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA,
	      double* X, int* INCX, double* BETA, double* Y, int* INCY);

  // Matrix-matrix multiply:
  // void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA,
  //	      double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);
}

using namespace booz_xform;

// Default constructor: set size to 1 x 1
Matrix::Matrix()
  : std::valarray<boozfloat>(1) // Call constructor of base class.
{
  nrows_ = 1;
  ncols_ = 1;
  len_ = 1;
}

Matrix::Matrix(index_type nrows_in, index_type ncols_in)
  : std::valarray<boozfloat>(nrows_in * ncols_in) // Call constructor of base class.
{
  nrows_ = nrows_in;
  ncols_ = ncols_in;
  len_ = nrows_ * ncols_;
}

void Matrix::resize(index_type nrows_in, index_type ncols_in, boozfloat v) {
  nrows_ = nrows_in;
  ncols_ = ncols_in;
  len_ = nrows_ * ncols_;
  std::valarray<boozfloat>::resize(nrows_ * ncols_, v);
}
	       
std::ostream& booz_xform::operator<< (std::ostream& os, Vector& v) {
  for (index_type j = 0; j < v.size(); j++) {
    if (j > 0) os << " ";
    os << v[j];
  }
  return os;
}

std::ostream& booz_xform::operator<< (std::ostream& os, Matrix& m) {
  for (index_type j = 0; j < m.nrows(); j++) {
    for (index_type k = 0; k < m.ncols(); k++) {
      if (k > 0) os << " ";
      os << m[j + m.nrows() * k];
    }
    os << std::endl;
  }
  return os;
}

boozfloat booz_xform::dot_product(Vector& u, Vector& v) {
  assert(u.size() == v.size());
  boozfloat sum = 0;
  for (int j = 0; j < u.size(); j++) {
    sum += u[j] * v[j];
  }
  return sum;
}


void booz_xform::matrix_vector_product(Matrix& m, Vector& v, Vector& result) {
  assert(m.ncols() == v.size());
  assert(m.nrows() == result.size());
  /*
  // Non-BLAS simple-minded method:
  boozfloat sum = 0;
  for (int j = 0; j < m.nrows(); j++) {
    sum = 0;
    for (int k = 0; k < v.size(); k++) sum += m(j, k) * v[k];
    result[j] = sum;
  }
  */

  char TRANS = 'N';
  int INC = 1;
  int nrows = m.nrows();
  int ncols = m.ncols();
  boozfloat ALPHA = 1.0;
  boozfloat BETA = 0.0;
  // dgemv would need to get replaced by sgemv if we switch to single precision
  dgemv_(&TRANS, &nrows, &ncols, &ALPHA, &m(0, 0), &nrows,
	      &v[0], &INC, &BETA, &result[0], &INC);
}

