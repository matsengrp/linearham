#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "subroutines.hpp"


TEST_CASE("colVecMatCwise", "[linalg]") {
  Eigen::VectorXd b(3);
  Eigen::MatrixXd A(3,4), B(3,4), correct_B(3,4);
  A << 1, 2.9, 3,  4,
       5, 6,   7,  8,
       9, 10,  11, 12;
  b << 0, 4, 1;
  correct_B <<  0,  0,  0,  0,
                20, 24, 28, 32,
                9, 10,  11, 12;

  colVecMatCwise(B, b, A);

  REQUIRE(B == correct_B);
}


TEST_CASE("rowVecMatCwise", "[linalg]") {
  Eigen::RowVectorXd b(4);
  Eigen::MatrixXd A(3,4), B(3,4), correct_B(3,4);
  A << 1, 2.9, 3,  4,
       5, 6,   7,  8,
       9, 10,  11, 12;
  b << 0, 4, 1, 10;
  correct_B <<  0, 11.6, 3, 40,
                0, 24,   7, 80,
                0, 40,  11, 120;

  rowVecMatCwise(B, b, A);

  REQUIRE(B == correct_B);
}


TEST_CASE("subProductMatrix", "[linalg]") {
  Eigen::MatrixXd A(3,3), correct_A(3,3);
  Eigen::VectorXd e(3);
  correct_A <<  -5,   -2, 2,
                -2.5, -1, 1,
                2.5,   1, 1;
  e << 2.5, -1, 2;
  A.setConstant(999);

  subProductMatrix(A, e);

  REQUIRE(A == correct_A);
}
