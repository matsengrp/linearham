#ifndef LINEARHAM_LINALG_
#define LINEARHAM_LINALG_

#include <Eigen/Dense>
#include <iostream>

/// @file linalg.hpp
/// @brief Linear algebra routines.

namespace linearham {


void ColVecMatCwise(const Eigen::Ref<const Eigen::VectorXd>& b,
                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                    Eigen::Ref<Eigen::MatrixXd> B);

void RowVecMatCwise(const Eigen::Ref<const Eigen::RowVectorXd>& b,
                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                    Eigen::Ref<Eigen::MatrixXd> B);

void SubProductMatrix(const Eigen::Ref<const Eigen::VectorXd>& e,
                      Eigen::Ref<Eigen::MatrixXd> A);
                      
void VectorByIndices(const Eigen::Ref<const Eigen::MatrixXd>& A,
                     const Eigen::Ref<const Eigen::VectorXi>& a,
                     Eigen::Ref<Eigen::VectorXd> b);

void BinaryMax(const Eigen::Ref<const Eigen::MatrixXd>& A,
               const Eigen::Ref<const Eigen::MatrixXd>& B,
               Eigen::Ref<Eigen::MatrixXd> C,
               Eigen::Ref<Eigen::MatrixXi> C_idx);
}

#endif  // LINEARHAM_LINALG_
