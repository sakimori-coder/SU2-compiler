#pragma once

# include <numeric>
#include <utility>
#include <vector>
#include <Eigen/Core>
#include <Eigen/LU>
#include "type.hpp"

namespace su2compiler 
{
namespace sdp
{


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
        for(int i = 0; i < M.size(); i++) {
            if(_M[i].cols() != _M[i].rows()) {
                throw std::domain_error("DiagBlockMatrix() : each block matrix must be square matrix");
            }
            BlockSizes[i] = _M[i].cols();
        }
        TotalSize = std::reduce(BlockSizes.begin(), BlockSizes.end());
    }

    DiagBlockMatrix(const std::vector<int>& _BlockSizes) {
        NumBlocks = _BlockSizes.size();
        BlockSizes = _BlockSizes;
        TotalSize = std::reduce(BlockSizes.begin(), BlockSizes.end());
        for(int i = 0; i < NumBlocks; i++) M[i] = MatrixXR::Zero(BlockSizes[i], BlockSizes[i]);
    }

    inline std::vector<int> structure() const noexcept { return BlockSizes; }
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
