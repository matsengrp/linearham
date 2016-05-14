#include <Eigen/Dense>

void ColVecMatCwise(
    Eigen::MatrixXd& B,
    Eigen::VectorXd& b,
    Eigen::MatrixXd& A);

void RowVecMatCwise(
    Eigen::MatrixXd& B,
    Eigen::RowVectorXd& b,
    Eigen::MatrixXd& A);

void SubProductMatrix(
    Eigen::MatrixXd& A,
    Eigen::VectorXd& e);

