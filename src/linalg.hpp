#include <Eigen/Dense>

void ColVecMatCwise(
    Eigen::Ref<Eigen::MatrixXd> B,
    const Eigen::Ref<const Eigen::VectorXd> b,
    const Eigen::Ref<const Eigen::MatrixXd> A);

void RowVecMatCwise(
    Eigen::Ref<Eigen::MatrixXd> B,
    const Eigen::Ref<const Eigen::RowVectorXd> b,
    const Eigen::Ref<const Eigen::MatrixXd> A);

void SubProductMatrix(
    Eigen::Ref<Eigen::MatrixXd> A,
    const Eigen::Ref<const Eigen::VectorXd> e);

