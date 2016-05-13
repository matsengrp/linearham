#include <Eigen/Dense>

void colVecMatCwise(
    Eigen::MatrixXd& B,
    Eigen::VectorXd& b,
    Eigen::MatrixXd& A);

void rowVecMatCwise(
    Eigen::MatrixXd& B,
    Eigen::RowVectorXd& b,
    Eigen::MatrixXd& A);

void subProductMatrix(
    Eigen::MatrixXd& A,
    Eigen::VectorXd& e);
