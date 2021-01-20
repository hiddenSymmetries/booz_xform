#ifndef BOOZ_VECTOR_MATRIX_H
#define BOOZ_VECTOR_MATRIX_H

#include <valarray>
#include <iostream>

namespace booz_xform {

  typedef double boozfloat;
  
  typedef std::valarray<boozfloat> Vector;
  
  class Matrix : public std::valarray<boozfloat> {

  private:
    size_t nrows_, ncols_, len_;

  public:
    //! Default constructor
    /**
     * Initialize a matrix with 1 row, 1 column, and a single entry 0.0.
     */
    Matrix();

    //! Constructor: Initialize a matrix and set all entries 0.0.
    /**
     * @param[in] nrows The number of rows desired.
     * @param[in] ncols The number of columns desired.
     */
    Matrix(size_t nrows, size_t ncols);

    //! Return the number of rows in the matrix.
    size_t nrows();

    //! Return the number of columns in the matrix.
    size_t ncols();

    //! Set the size of the matrix, allocating memory for it.
    /**
     * @param[in] nrows The number of rows desired.
     * @param[in] ncols The number of columns desired.
     * @param[in] value All the elements of the matrix will be set to this value.
     */
    void resize(size_t nrows, size_t ncols, boozfloat value);
    
    // For info about matrix indexing:
    // https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
    //! Access a specified element of the matrix
    /**
     * @param[in] row The 0-based index of the row to access.
     * @param[in] col The 0-based index of the column to access. 
     */
    boozfloat& operator()(size_t row, size_t column);
    boozfloat  operator()(size_t, size_t) const;

    Matrix& operator=(const boozfloat);
    Matrix& operator=(const Matrix&);
  };

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
  inline size_t Matrix::nrows() {
    return nrows_;
  }

  inline size_t Matrix::ncols() {
    return ncols_;
  }

  inline boozfloat& Matrix::operator()(size_t m, size_t n) {
    return (*this)[m + nrows_ * n];
  }

  inline boozfloat Matrix::operator()(size_t m, size_t n) const {
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
