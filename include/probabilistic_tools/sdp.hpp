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

using DenseMat  = MatrixXR;
using SparseMat = Eigen::SparseMatrix<Real>;
using DenseVec  = VectorXR;
using SparseVec = Eigen::SparseVector<Real>;
using MatVecVar = std::variant<DenseMat,
                               SparseMat,
                               DenseVec,
                               SparseVec>;




class DiagBlockMatrix
{
private:
    std::vector<MatVecVar> Blocks;
    std::vector<int>       BlockSizes;
    UINT                   NumBlocks = 0;
    UINT                   TotalSize = 0;

    

public:
    DiagBlockMatrix(const std::vector<MatVecVar>& _Blocks) {
        Blocks = _Blocks;
        NumBlocks = Blocks.size();
        BlockSizes.resize(NumBlocks);

        for(int i = 0; i < NumBlocks; i++) {
            std::visit([&](const auto& blk) {
                using T = std::decay_t<decltype(blk)>;

                if constexpr (std::is_same_v<T, DenseMat> || std::is_same_v<T, SparseMat>) {
                    if(blk.rows() != blk.cols()) {
                        throw std::domain_error("DiagBlockMatrix() : block matrix must be square matrix");
                    }
                    BlockSizes[i] = blk.rows();
                }
                else if constexpr (std::is_same_v<T, DenseVec> || std::is_same_v<T, SparseVec>) {
                    BlockSizes[i] = -blk.size();
                }
            }, Blocks[i]);
        }

        TotalSize = 0;
        for(int size : BlockSizes) TotalSize += std::abs(size);
    }


    inline std::vector<MatVecVar> get_Blocks() const noexcept { return Blocks; }
    inline std::vector<int> get_BlockSizes() const noexcept { return BlockSizes; }
    inline UINT get_TotalSize() const noexcept { return TotalSize; }
    MatrixXR to_DenseMatrix() const;
    DiagBlockMatrix transpose() const;
    DiagBlockMatrix inverse() const;
    Real maxAbsCoeff() const;
    DiagBlockMatrix LLT() const;
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
