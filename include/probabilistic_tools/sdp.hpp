#pragma once

# include <numeric>
#include <utility>
#include <vector>
#include <variant>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include "type.hpp"

namespace su2compiler 
{
namespace sdp
{

using MatrixVar = std::variant<MatrixXR, Eigen::SparseMatrix<Real>>;

class DiagBlockMatrix
{
private:
    int TotalSize;
    int NumBlocks;
    std::vector<int> BlockSizes;
    std::vector<MatrixXR> M;
    

public:
    DiagBlockMatrix(const std::vector<MatrixXR>& _M) {
        M = _M;
        NumBlocks = M.size();
        BlockSizes.resize(NumBlocks);
        for(int i = 0; i < NumBlocks; i++) {
            if(_M[i].cols() == 1) {
                BlockSizes[i] = -_M[i].rows();
            } else if(_M[i].rows() == _M[i].cols()) {
                BlockSizes[i] = _M[i].rows();
            } else {
                throw std::domain_error("DiagBlockMatrix() : each block matrix must be square matrix");
            }
        }
        TotalSize = 0;
        for(int size : BlockSizes) TotalSize += std::abs(size);
    }

    DiagBlockMatrix(const std::vector<int>& _BlockSizes) {
        NumBlocks = _BlockSizes.size();
        BlockSizes = _BlockSizes;
        TotalSize = 0;
        for(int size : BlockSizes) TotalSize += std::abs(size);
        for(int i = 0; i < NumBlocks; i++) {
            if(BlockSizes[i] > 0) M[i] = MatrixXR::Zero(BlockSizes[i], BlockSizes[i]);
            else                  M[i] = VectorXR::Zero(-BlockSizes[i]);
        }
    }

    inline std::vector<int> structure() const noexcept { return BlockSizes; }
    inline int get_TotalSize() const noexcept { return TotalSize; }
    MatrixXR to_DenseMatrix() const noexcept;

    DiagBlockMatrix transpose() const noexcept;
    DiagBlockMatrix inverse() const;
    Real maxAbsCoeff() const noexcept;
    DiagBlockMatrix cholesky_decomposition() const;
    std::vector<Real> compute_eigenvalues() const;


    DiagBlockMatrix& operator+=(const DiagBlockMatrix& r);
    DiagBlockMatrix& operator-=(const DiagBlockMatrix& r);
    DiagBlockMatrix& operator*=(const DiagBlockMatrix& r);
    DiagBlockMatrix& operator*=(const Real& r) noexcept;

    DiagBlockMatrix operator-() const noexcept;

    friend std::ostream& operator<<(std::ostream& os, const DiagBlockMatrix& x);
    friend Real HSinner(const DiagBlockMatrix& A, const DiagBlockMatrix& B);
    friend bool isSymmetric(const DiagBlockMatrix& A, Real tol);
};

inline DiagBlockMatrix operator+(DiagBlockMatrix lhs, const DiagBlockMatrix& rhs) { return lhs += rhs; }
inline DiagBlockMatrix operator-(DiagBlockMatrix lhs, const DiagBlockMatrix& rhs) { return lhs -= rhs; }
inline DiagBlockMatrix operator*(DiagBlockMatrix lhs, const DiagBlockMatrix& rhs) { return lhs *= rhs; }
inline DiagBlockMatrix operator*(DiagBlockMatrix lhs, const Real& rhs) { return lhs *= rhs; }
inline DiagBlockMatrix operator*(const Real lhs, DiagBlockMatrix rhs) { return rhs *= lhs; }

DiagBlockMatrix BlockIdentity(const std::vector<int>& structure);
Real HSinner(const DiagBlockMatrix& A, const DiagBlockMatrix& B);    

bool isSymmetric(const MatrixXR& A, Real tol);
bool isSymmetric(const DiagBlockMatrix& A, Real tol);


VectorXR SDP(
    const VectorXR& c,
    const std::vector<DiagBlockMatrix>& F, 
    Real lambda,
    Real betaStar,
    Real betaBar,
    Real gamma,
    Real epsilon1,
    Real epsilon2,
    int maxITERATION,
    bool OUTPUT_HISTORY
);

}
}
