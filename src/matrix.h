/*
 * Matrix
 * 
 * This library implements the matrix
 * data type and math operations with it
 * 
 * Copyright:
 * Eldar dougiela Bagirov
 * 
 */

#ifndef _MATRIX_LIB
#define _MATRIX_LIB

#define _equality_eps 1e-7

class Matrix {
 private:
  int rows_;
  int cols_;
  double* matrix_;

  static double* _allocateMatrix(const int r, const int c);
  static void _freeMatrix(double*& m);
  static void _copyMatrix(Matrix& dst, const Matrix& src);

  double _getMatrixElementForDeterminant(const int r, const int c) const;
  void _swapRows(const int a, const int b);

 public:
  Matrix();
  Matrix(const int rows, const int cols);
  Matrix(const Matrix& other);
  Matrix(Matrix&& other);
  ~Matrix();

  int getCols() const;
  int getRows() const;
  double* getMatrixPtr();
  void setCols(const int newCols);
  void setRows(const int newRows);
  void setMatrixPtr(double* newBuff);

  bool EqMatrix(const Matrix&) const;
  void SumMatrix(const Matrix&);
  void SubMatrix(const Matrix&);
  void MulNumber(const double);
  void MulMatrix(const Matrix&);

  Matrix Transpose() const;
  Matrix CalcComplements() const;
  double getMinor(const int r, const int c) const;
  Matrix exclude(const int r, const int c) const;
  double Determinant() const;
  double modifyingDeterminant();
  Matrix InverseMatrix() const;

  Matrix operator+(const Matrix&);
  Matrix operator-(const Matrix&);
  Matrix operator*(const Matrix&);
  Matrix operator*(const double);
  bool operator==(const Matrix&) const;
  bool operator!=(const Matrix&) const;
  Matrix& operator=(const Matrix&);
  Matrix& operator=(Matrix&&);
  void operator+=(const Matrix&);
  void operator-=(const Matrix&);
  void operator*=(const Matrix&);
  void operator*=(const double);
  double& operator()(const int r, const int c);
  const double& operator()(const int r, const int c) const;
};

#endif
