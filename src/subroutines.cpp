#include <Eigen/Dense>
#include <iostream>
using namespace std;


void rowVecMatCwise(
    Eigen::MatrixXd& A,
    Eigen::RowVectorXd& b,
    Eigen::MatrixXd& B) {
  for(int i=0; i < A.rows(); i++) {
    A.row(i) = b.cwiseProduct(B.row(i));
  }
}

/*
int main()
{
  Eigen::MatrixXd A(3,4);
  Eigen::RowVectorXd b(4);
  Eigen::VectorXd c(4);
  Eigen::MatrixXd B(3,4);
  B <<  1, 2.9, 3, 4,
        5, 6, 7, 8,
        9,10,11,12;
  b << 0, 4, 1, 10;

  test(A, b, B);
  cout << "orig\n" << B << endl;
  cout << "final\n" << A << endl;

}
*/
