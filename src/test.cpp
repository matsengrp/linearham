#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "subroutines.hpp"


TEST_CASE("Componentwise", "[componentwise]") {
  Eigen::RowVectorXd b(4);
  Eigen::MatrixXd A(3,4), B(3,4), correct_A(3,4);
  B <<  1, 2.9, 3, 4,
        5, 6, 7, 8,
        9,10,11,12;
  b << 0, 4, 1, 10;
  correct_A <<  0, 11.6, 3, 40,
                0, 24,   7, 80,
                0, 40,  11, 120;

  rowVecMatCwise(A, b, B);

  REQUIRE(A == correct_A);
}
