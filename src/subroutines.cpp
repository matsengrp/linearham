#include <Eigen/Dense>
#include <iostream>
using namespace std;


void rowVecMatCwise(
    Eigen::MatrixXd& B,
    Eigen::RowVectorXd& b,
    Eigen::MatrixXd& A) {
  for(int i=0; i < B.rows(); i++) {
    B.row(i) = b.cwiseProduct(A.row(i));
  }
}
