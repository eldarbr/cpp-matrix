#include <gtest/gtest.h>

#include <cmath>

#include "../matrix.h"

TEST(complements, incorrect_matrix) {
  Matrix m(10, 20);
  EXPECT_THROW(m.CalcComplements(), std::invalid_argument);
}

TEST(complements, empty_matrix) {
  Matrix m;
  EXPECT_THROW(m.CalcComplements(), std::invalid_argument);
}

TEST(complements, successful_1) {
  Matrix m(3, 3);
  m(0, 0) = 1;
  m(0, 1) = 2;
  m(0, 2) = 3;
  m(1, 0) = 0;
  m(1, 1) = 4;
  m(1, 2) = 2;
  m(2, 0) = 5;
  m(2, 1) = 2;
  m(2, 2) = 1;

  Matrix expected(3, 3);
  expected(0, 0) = 0;
  expected(0, 1) = 10;
  expected(0, 2) = -20;
  expected(1, 0) = 4;
  expected(1, 1) = -14;
  expected(1, 2) = 8;
  expected(2, 0) = -8;
  expected(2, 1) = -2;
  expected(2, 2) = 4;

  Matrix actual;
  EXPECT_NO_THROW(actual = m.CalcComplements());
  EXPECT_EQ(actual, expected);
}

TEST(complements, successful_2) {
  Matrix m(3, 3);
  m(0, 0) = 1;
  m(0, 1) = 2;
  m(0, 2) = 3;
  m(1, 1) = 4;
  m(1, 2) = 2;
  m(2, 0) = 5;
  m(2, 1) = 2;
  m(2, 2) = 1;

  Matrix expected(3, 3);
  expected(0, 1) = 10;
  expected(0, 2) = -20;
  expected(1, 0) = 4;
  expected(1, 1) = -14;
  expected(1, 2) = 8;
  expected(2, 0) = -8;
  expected(2, 1) = -2;
  expected(2, 2) = 4;

  Matrix actual;
  EXPECT_NO_THROW(actual = m.CalcComplements());
  EXPECT_EQ(actual, expected);
}

TEST(determinant, non_initialized) {
  Matrix m;
  EXPECT_THROW(m.CalcComplements(), std::invalid_argument);
}

TEST(determinant, non_square) {
  Matrix m(1, 2);
  EXPECT_THROW(m.CalcComplements(), std::invalid_argument);
}

TEST(determinant, successful_1) {
  Matrix m(3, 3);
  m(0, 0) = 1;
  m(0, 1) = 2;
  m(0, 2) = 3;
  m(1, 0) = 0;
  m(1, 1) = 4;
  m(1, 2) = 2;
  m(2, 0) = 5;
  m(2, 1) = 2;
  m(2, 2) = 1;

  const double expected = -40.;

  double actual;
  EXPECT_NO_THROW(actual = m.Determinant());
  EXPECT_FLOAT_EQ(actual, expected);
}

TEST(determinant, successful_2) {
  Matrix m(4, 4);
  m(0, 0) = 5;
  m(0, 1) = 12;
  m(0, 2) = 7;
  m(0, 3) = 3;
  m(1, 0) = 9;
  m(1, 1) = 15;
  m(1, 2) = 8;
  m(1, 3) = 6;
  m(2, 0) = 4;
  m(2, 1) = 1;
  m(2, 2) = 11;
  m(2, 3) = 14;
  m(3, 0) = 2;
  m(3, 1) = 10;
  m(3, 2) = 13;
  m(3, 3) = 17;

  const double expected = -3771.;

  double actual;
  EXPECT_NO_THROW(actual = m.Determinant());
  EXPECT_FLOAT_EQ(actual, expected);
}

TEST(eq, nonexistent_matrices) {
  Matrix m1;
  Matrix m2;
  EXPECT_EQ(m1, m2);
}

TEST(eq, equal_empty_matrices) {
  Matrix m1(10, 20);
  Matrix m2(10, 20);
  EXPECT_EQ(m1, m2);
}

TEST(eq, not_equal_empty_matrices) {
  Matrix m1(10, 20);
  Matrix m2(20, 20);
  EXPECT_NE(m1, m2);
}

TEST(eq, not_equal_not_empty_matrices) {
  Matrix m1(20, 20);
  m1(0, 5) = 12.;
  m1(11, 8) = -2.;
  Matrix m2(20, 20);
  m2(3, 5) = 12.;
  m2(11, 8) = -2.;
  EXPECT_NE(m1, m2);
}

TEST(eq, equal_not_empty_matrices) {
  Matrix m1(20, 20);
  m1(3, 5) = 12.;
  m1(11, 8) = -2.;
  Matrix m2(20, 20);
  m2(3, 5) = 12.;
  m2(11, 8) = -2.;
  EXPECT_EQ(m1, m2);
}

TEST(inverse, incorrect) {
  Matrix m;
  EXPECT_THROW(m.InverseMatrix(), std::invalid_argument);
}

TEST(inverse, non_square) {
  Matrix m(1, 2);
  EXPECT_THROW(m.InverseMatrix(), std::invalid_argument);
}

TEST(inverse, successful_1) {
  Matrix m1(3, 3);
  m1(0, 0) = 1;
  m1(0, 1) = 2;
  m1(0, 2) = 3;
  m1(1, 0) = 0;
  m1(1, 1) = 4;
  m1(1, 2) = 2;
  m1(2, 0) = 5;
  m1(2, 1) = 2;
  m1(2, 2) = 1;

  Matrix expected(3, 3);
  expected(0, 0) = 0;
  expected(0, 1) = -.1;
  expected(0, 2) = .2;
  expected(1, 0) = -.25;
  expected(1, 1) = .35;
  expected(1, 2) = .05;
  expected(2, 0) = .5;
  expected(2, 1) = -.2;
  expected(2, 2) = -.1;

  Matrix actual;
  EXPECT_NO_THROW(actual = m1.InverseMatrix());
  EXPECT_EQ(actual, expected);
}

TEST(inverse, successful_2) {
  Matrix m1(4, 4);
  m1(0, 0) = 5;
  m1(0, 1) = 12;
  m1(0, 2) = 7;
  m1(0, 3) = 3;
  m1(1, 0) = 9;
  m1(1, 1) = 15;
  m1(1, 2) = 8;
  m1(1, 3) = 6;
  m1(2, 0) = 4;
  m1(2, 1) = 1;
  m1(2, 2) = 11;
  m1(2, 3) = 14;
  m1(3, 0) = 2;
  m1(3, 1) = 10;
  m1(3, 2) = 13;
  m1(3, 3) = 17;

  Matrix expected(4, 4);
  expected(0, 0) = -0.126492;
  expected(0, 1) = 0.167064;
  expected(0, 2) = 0.0859189;
  expected(0, 3) = -0.107399;
  expected(1, 0) = -0.0251923;
  expected(1, 1) = 0.043755;
  expected(1, 2) = -0.104482;
  expected(1, 3) = 0.0750464;
  expected(2, 0) = 0.392204;
  expected(2, 1) = -0.260143;
  expected(2, 2) = 0.142403;
  expected(2, 3) = -0.0946698;
  expected(3, 0) = -0.27022;
  expected(3, 1) = 0.15354;
  expected(3, 2) = -0.0575444;
  expected(3, 3) = 0.0997083;

  Matrix actual;
  EXPECT_NO_THROW(actual = m1.InverseMatrix());
  for (int i = 0; i != actual.getRows(); ++i) {
    for (int j = 0; j != actual.getCols(); ++j) {
      EXPECT_LE(fabs(expected(i, j) - actual(i, j)), 1e-6);
    }
  }
}

TEST(mult_num, uninitialized) {
  Matrix m;
  EXPECT_THROW(m *= 1, std::invalid_argument);
}

TEST(mult_num, successful) {
  const double multiplier = 12;
  Matrix m1(10, 2);
  m1(5, 0) = 10;
  m1(1, 1) = -15;
  Matrix expected(10, 2);
  expected(5, 0) = 120;
  expected(1, 1) = -180;

  Matrix actual;
  EXPECT_NO_THROW(actual = m1 * multiplier);

  EXPECT_EQ(actual, expected);
}

TEST(sum, successful) {
  Matrix m1(10, 2);
  Matrix m2(10, 2);
  m1(5, 0) = 10;
  m2(1, 1) = -15;
  Matrix expected(10, 2);
  expected(5, 0) = 10;
  expected(1, 1) = -30;

  Matrix actual;
  EXPECT_NO_THROW(actual = m1 + m2);
  EXPECT_NO_THROW(actual += m2);

  EXPECT_EQ(actual, expected);
}

TEST(sub, successful) {
  Matrix m1(10, 2);
  Matrix m2(10, 2);
  m1(5, 0) = 10;
  m2(1, 1) = -15;
  Matrix expected(10, 2);
  expected(5, 0) = 10;
  expected(1, 1) = 30;

  Matrix actual;
  EXPECT_NO_THROW(actual = m1 - m2);
  EXPECT_NO_THROW(actual -= m2);

  EXPECT_EQ(actual, expected);
}

TEST(move, successful) {
  Matrix m1(10, 2);
  Matrix m2(10, 2);
  m1(5, 0) = 10;
  m2(1, 1) = -15;

  Matrix r1(std::move(m1));
  EXPECT_EQ(m1.getMatrixPtr(), nullptr);
  EXPECT_EQ(r1(5, 0), 10);

  EXPECT_NO_THROW(r1 = std::move(m2));
  EXPECT_EQ(m2.getMatrixPtr(), nullptr);
  EXPECT_EQ(r1(5, 0), 0);
  EXPECT_EQ(r1(1, 1), -15);
}

TEST(mult, successful) {
  Matrix m1(3, 2);
  m1(0, 0) = 1, m1(0, 1) = 4;
  m1(1, 0) = 2, m1(1, 1) = 5;
  m1(2, 0) = 3, m1(2, 1) = 6;

  Matrix m2(2, 3);
  m2(0, 0) = 1, m2(0, 1) = -1, m2(0, 2) = 1;
  m2(1, 0) = 2, m2(1, 1) = 3, m2(1, 2) = 4;

  Matrix expected(3, 3);
  expected(0, 0) = 9;
  expected(0, 1) = 11;
  expected(0, 2) = 17;
  expected(1, 0) = 12;
  expected(1, 1) = 13;
  expected(1, 2) = 22;
  expected(2, 0) = 15;
  expected(2, 1) = 15;
  expected(2, 2) = 27;

  Matrix result;
  EXPECT_NO_THROW(result = m1 * m2);
  EXPECT_EQ(result, expected);

  EXPECT_THROW(result *= m2, std::invalid_argument);
}

int main(int argc, char* argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
