#ifndef BOOZ_VECTOR_MATRIX_H
#define BOOZ_VECTOR_MATRIX_H

#include <valarray>
#include <iostream>

namespace booz_xform {

  typedef double boozfloat;
  
  typedef std::valarray<boozfloat> Vector;
  
  // typedef unsigned index_type;
  typedef int index_type;

  class Matrix : public std::valarray<boozfloat> {

  private:
    index_type nrows_, ncols_, len_;

  public:
    Matrix();
    Matrix(index_type, index_type);
    index_type nrows();
    index_type ncols();
    void resize(index_type, index_type, boozfloat);
    // For info about matrix indexing:
    // https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
    boozfloat& operator()(index_type, index_type);
    boozfloat  operator()(index_type, index_type) const;
    Matrix& operator=(const boozfloat);
    Matrix& operator=(const Matrix&);
  };

  void matrix_vector_product(Matrix&, Vector&, Vector&);
  boozfloat dot_product(Vector&, Vector&);
  
  // These operators should be in the namespace:
  // https://stackoverflow.com/questions/3891402/operator-overloading-and-namespaces
  // https://stackoverflow.com/questions/3623631/where-should-non-member-operator-overloads-be-placed
  std::ostream& operator<<(std::ostream&, Vector&);
  std::ostream& operator<<(std::ostream&, Matrix&);
  
  // Multiplication of an int with std::valarray<boozfloat> is not included in some compilers:
  inline Vector operator*(int j, Vector& v) {
    //return boozfloat(j) * v;
    return std::operator*(boozfloat(j), v);
  }

  // gcc complains about ambiguity if we do not have this next little function:
  inline Vector operator*(boozfloat j, Vector& v) {
    return std::operator*(j, v);
  }
  
  /*
  inline Vector operator*(float j, Vector& v) {
    return std::operator*(boozfloat(j), v);
  }
  
  inline Vector operator*(double j, Vector& v) {
    return std::operator*(boozfloat(j), v);
  }
  */

  // inline functions must be included in every file that uses them,
  // so these functions should go in the header file.
  inline index_type Matrix::nrows() {
    return nrows_;
  }

  inline index_type Matrix::ncols() {
    return ncols_;
  }

  inline boozfloat& Matrix::operator()(index_type m, index_type n) {
    return (*this)[m + nrows_ * n];
  }

  inline boozfloat Matrix::operator()(index_type m, index_type n) const {
    return (*this)[m + nrows_ * n];
  }

  inline Matrix& Matrix::operator=(const boozfloat v) {
    // Delegate to parent class:
    std::valarray<boozfloat>::operator=(v);
    return *this;
  }

  inline Matrix& Matrix::operator=(const Matrix &v) {
    // Delegate to parent class:
    std::valarray<boozfloat>::operator=(v);
    return *this;
  }


} // namespace booz_xform


#endif
