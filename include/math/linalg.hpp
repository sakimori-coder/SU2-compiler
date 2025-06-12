#pragma once

#include <Eigen/Core>

#include "core/type.hpp"


namespace su2compiler::math::linalg {

MatrixR GSO(const MatrixR& B);

// compute lower triangular matrix L such that A = L * L^T
MatrixR Cholesky(const MatrixR& A);


VectorR SolveSystemSPD(
        const MatrixR& A,
        const VectorR& b
);


MatrixR InverseSPD(const MatrixR& A);

}