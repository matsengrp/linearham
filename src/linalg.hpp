#ifndef LINEARHAM_LINALG_
#define LINEARHAM_LINALG_

#include <Eigen/Dense>

/// @file linalg.hpp
/// @brief Linear algebra routines.

namespace linearham {

  const double SCALE_FACTOR = pow(2, 256);
  const double SCALE_THRESHOLD = (1.0 / SCALE_FACTOR);
  const double LOG_SCALE_FACTOR = log(SCALE_FACTOR);

void ColVecMatCwise(const Eigen::Ref<const Eigen::VectorXd>& b,
                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                    Eigen::Ref<Eigen::MatrixXd> B);

void RowVecMatCwise(const Eigen::Ref<const Eigen::RowVectorXd>& b,
                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                    Eigen::Ref<Eigen::MatrixXd> B);
int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m);
int SubProductMatrix(const Eigen::Ref<const Eigen::VectorXd>& e,
                      Eigen::Ref<Eigen::MatrixXd> A);

void VectorByIndices(const Eigen::Ref<const Eigen::MatrixXd>& A,
                     const Eigen::Ref<const Eigen::VectorXi>& a,
                     Eigen::Ref<Eigen::VectorXd> b);

void BinaryMax(const Eigen::Ref<const Eigen::MatrixXd>& A,
               const Eigen::Ref<const Eigen::MatrixXd>& B,
               Eigen::Ref<Eigen::MatrixXd> C,
               Eigen::Ref<Eigen::MatrixXi> C_idx);


}  // namespace linearham

#endif  // LINEARHAM_LINALG_
