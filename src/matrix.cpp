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

#include "matrix.h"

#include <cmath>
#include <cstring>
#include <stdexcept>
#include <utility>

//////////////////////////////////////////////////
//                    Matrix                    //
//            begin of implementation           //

Matrix::Matrix() {
  cols_ = 0;
  rows_ = 0;
  matrix_ = nullptr;
}

Matrix::Matrix(const int r, const int c) : rows_(r), cols_(c) {
  matrix_ = _allocateMatrix(rows_, cols_);
}

Matrix::~Matrix() { _freeMatrix(matrix_); }

Matrix::Matrix(const Matrix& src) : matrix_(nullptr) { *this = src; }

Matrix::Matrix(Matrix&& src) : matrix_(nullptr) { *this = std::move(src); }

int Matrix::getCols() const { return cols_; }

int Matrix::getRows() const { return rows_; }

double* Matrix::getMatrixPtr() { return matrix_; }

void Matrix::setCols(const int new_c) {
  if (new_c) {
    Matrix new_m(rows_, new_c);
    _copyMatrix(new_m, *this);
    *this = std::move(new_m);
  } else {
    _freeMatrix(matrix_);
    matrix_ = nullptr;
    rows_ = 0;
    cols_ = 0;
  }
}

void Matrix::setRows(const int new_r) {
  if (new_r) {
    Matrix new_m(new_r, cols_);
    _copyMatrix(new_m, *this);
    *this = std::move(new_m);
  } else {
    _freeMatrix(matrix_);
    matrix_ = nullptr;
    rows_ = 0;
    cols_ = 0;
  }
}

void Matrix::setMatrixPtr(double* matrix) { this->matrix_ = matrix; }

bool Matrix::EqMatrix(const Matrix& other) const {
  bool equal = true;
  if (rows_ != other.getRows() || cols_ != other.getCols()) {
    equal = false;
  }
  for (int i = 0; equal && i < rows_; ++i) {
    for (int j = 0; equal && j < cols_; ++j) {
      const double delta = fabs((*this)(i, j) - other(i, j));
      if (delta > _equality_eps) {
        equal = false;
      }
    }
  }
  return equal;
}

void Matrix::SumMatrix(const Matrix& other) {
  if (rows_ != other.getRows() || cols_ != other.getCols()) {
    throw std::invalid_argument("wrong matrix sizes");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      (*this)(i, j) += other(i, j);
    }
  }
}

void Matrix::SubMatrix(const Matrix& other) {
  if (rows_ != other.getRows() || cols_ != other.getCols()) {
    throw std::invalid_argument("wrong matrix sizes");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      (*this)(i, j) -= other(i, j);
    }
  }
}

void Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      (*this)(i, j) *= num;
    }
  }
}

void Matrix::MulMatrix(const Matrix& other) {
  if (cols_ != other.getRows()) {
    throw std::invalid_argument("multiplication-incompatible matrix sizes");
  }
  Matrix res(rows_, other.getCols());
  for (int i = 0; i < res.getRows(); ++i) {
    for (int j = 0; j < res.getCols(); ++j) {
      for (int k = 0; k < cols_; ++k) {
        res(i, j) += (*this)(i, k) * other(k, j);
      }
    }
  }
  *this = std::move(res);
}

Matrix Matrix::Transpose() const {
  Matrix res(cols_, rows_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      res(j, i) = (*this)(i, j);
    }
  }
  return res;
}

Matrix Matrix::CalcComplements() const {
  if (cols_ != rows_) {
    throw std::invalid_argument("the matrix should be square");
  }
  Matrix result(rows_, cols_);
  const int side = rows_;
  for (int i = 0; i < side; ++i) {
    for (int j = 0; j < side; ++j) {
      const int sign = 1 - (((i + j) & 1) << 1);
      double minor = getMinor(i, j);
      result(i, j) = sign * minor;
    }
  }
  return result;
}

double Matrix::getMinor(const int excluded_row, const int excluded_col) const {
  if (rows_ != cols_) {
    throw std::invalid_argument(
        "the matrix should be square for minor to be calculated");
  }
  if (cols_ < 2) {
    return 1.;
  }
  Matrix tmp = exclude(excluded_row, excluded_col);
  return tmp.modifyingDeterminant();
}

Matrix Matrix::exclude(const int r, const int c) const {
  if (cols_ < 2 || rows_ < 2) {
    throw std::invalid_argument("too small to exclude");
  }
  Matrix result(rows_ - 1, cols_ - 1);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (i != r && j != c) {
        const int dst_i = i - (i > r);
        const int dst_j = j - (j > c);
        result(dst_i, dst_j) = (*this)(i, j);
      }
    }
  }
  return result;
}

double Matrix::Determinant() const {
  if (cols_ != rows_) {
    throw std::invalid_argument("the matrix should be square");
  }
  Matrix copy(*this);
  const double determinant = copy.modifyingDeterminant();
  return determinant;
}

double Matrix::modifyingDeterminant() {
  if (cols_ != rows_) {
    throw std::invalid_argument("the matrix should be square");
  }
  double sign = 1;
  const int side = rows_;
  int zero_divisor = 0;
  for (int k = 0; !zero_divisor && k < side - 1; ++k) {
    if (fabs(_getMatrixElementForDeterminant(k, k)) < 1e-9) {
      int next_swap = side;
      for (int m = k + 1; next_swap == side && m < side; ++m) {
        if (fabs(_getMatrixElementForDeterminant(m, k)) >= 1e-9) {
          next_swap = m;
        }
      }
      if (next_swap == side) {
        zero_divisor = 1;
      } else {
        _swapRows(next_swap, k);
        sign = -sign;
      }
    }
    for (int i = k + 1; !zero_divisor && i < side; ++i) {
      for (int j = k + 1; !zero_divisor && j < side; ++j) {
        const double numerator = (_getMatrixElementForDeterminant(i, j) *
                                  _getMatrixElementForDeterminant(k, k)) -
                                 (_getMatrixElementForDeterminant(i, k) *
                                  _getMatrixElementForDeterminant(k, j));
        const double divisor = _getMatrixElementForDeterminant(k - 1, k - 1);
        if (divisor != 0) {
          (*this)(i, j) = numerator / divisor;
        } else {
          zero_divisor = 1;
        }
      }
    }
  }
  if (zero_divisor) {
    (*this)(side - 1, side - 1) = 0.;
  }
  const double result = sign * (*this)(side - 1, side - 1);
  return result;
}

double Matrix::_getMatrixElementForDeterminant(const int r, const int c) const {
  double element = 0.;
  if (r < 0 || c < 0) {
    element = 1.;
  } else if (r < rows_ && c < cols_) {
    element = (*this)(r, c);
  }
  return element;
}

void Matrix::_swapRows(const int a, const int b) {
  for (int i = 0; i < cols_; ++i) {
    const double tmp = (*this)(a, i);
    (*this)(a, i) = (*this)(b, i);
    (*this)(b, i) = tmp;
  }
}

Matrix Matrix::InverseMatrix() const {
  const double determinant = Determinant();
  if (fabs(determinant - 0.) < 1e-6) {
    throw std::invalid_argument(
        "inverse matrix is only avaliable for a nonsingular matrix");
  }
  Matrix transposedComplementsMatrix = CalcComplements().Transpose();
  transposedComplementsMatrix.MulNumber(1. / determinant);
  return transposedComplementsMatrix;
}

Matrix Matrix::operator+(const Matrix& other) {
  Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

Matrix Matrix::operator-(const Matrix& other) {
  Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

Matrix Matrix::operator*(const Matrix& other) {
  Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

Matrix Matrix::operator*(const double num) {
  Matrix result(*this);
  result.MulNumber(num);
  return result;
}

bool Matrix::operator==(const Matrix& other) const { return EqMatrix(other); }

bool Matrix::operator!=(const Matrix& other) const { return !EqMatrix(other); }

Matrix& Matrix::operator=(const Matrix& other) {
  if (this != &other) {
    _freeMatrix(matrix_);
    rows_ = other.getRows();
    cols_ = other.getCols();
    if (rows_ && cols_) {
      matrix_ = _allocateMatrix(rows_, cols_);
      _copyMatrix(*this, other);
    } else {
      matrix_ = nullptr;
    }
  }
  return *this;
}

Matrix& Matrix::operator=(Matrix&& src) {
  if (this != &src) {
    _freeMatrix(matrix_);
    rows_ = src.getRows();
    cols_ = src.getCols();
    matrix_ = src.getMatrixPtr();
    src.setMatrixPtr(nullptr);
    src.setCols(0);
    src.setRows(0);
  }
  return *this;
}

void Matrix::operator+=(const Matrix& other) {
  if (nullptr == matrix_) {
    throw std::invalid_argument("the matrix is uninitialized");
  }
  SumMatrix(other);
}

void Matrix::operator-=(const Matrix& other) {
  if (nullptr == matrix_) {
    throw std::invalid_argument("the matrix is uninitialized");
  }
  SubMatrix(other);
}

void Matrix::operator*=(const Matrix& other) {
  if (nullptr == matrix_) {
    throw std::invalid_argument("the matrix is uninitialized");
  }
  MulMatrix(other);
}

void Matrix::operator*=(const double num) {
  if (nullptr == matrix_) {
    throw std::invalid_argument("the matrix is uninitialized");
  }
  MulNumber(num);
}

double& Matrix::operator()(const int r, const int c) {
  if (r < 0 || c < 0 || r >= rows_ || c >= cols_) {
    throw std::invalid_argument("rows or cols not in range");
  }
  return matrix_[r * cols_ + c];
}

const double& Matrix::operator()(const int r, const int c) const {
  return (const double&)((*((Matrix*)this))(r, c));
}

double* Matrix::_allocateMatrix(const int r, const int c) {
  if (r < 1 || c < 1) {
    throw std::invalid_argument("rows or cols not in range");
  }
  double* m = new double[r * c];
  memset(m, 0, r * c * sizeof(double));
  return m;
}

void Matrix::_freeMatrix(double*& m) {
  if (nullptr != m) {
    delete[] m;
  }
}

void Matrix::_copyMatrix(Matrix& dst, const Matrix& src) {
  const int minRows =
      dst.getRows() < src.getRows() ? dst.getRows() : src.getRows();
  const int minCols =
      dst.getCols() < src.getCols() ? dst.getCols() : src.getCols();
  for (int i = 0; i < minRows; ++i) {
    for (int j = 0; j < minCols; ++j) {
      dst(i, j) = src(i, j);
    }
  }
}

//                    Matrix                    //
//             end of implementation            //
//////////////////////////////////////////////////
