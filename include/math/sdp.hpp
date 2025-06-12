#pragma once

#include <cmath>
#include <utility>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "core/type.hpp"

namespace su2compiler::math::sdp
{

using SparseMatrixR = Eigen::SparseMatrix<Real>;

class DenseDiagBlock
{
public:
    int NumBlocks = 0;
    std::vector<MatrixR> Blocks;
    std::vector<int> BlockSizes;   // If BlockSizes[i] < 0, i-th block is diagonal matrix
    

    DenseDiagBlock(std::vector<MatrixR>& _Blocks) {
        Blocks = _Blocks;
        NumBlocks = Blocks.size();
        for(auto blk : Blocks) {
            if(blk.cols() != 1) {
                assert(blk.rows() == blk.cols());
                BlockSizes.push_back( blk.rows());
            } else {
                BlockSizes.push_back(-blk.rows());
            }
        }
    }

    DenseDiagBlock(std::vector<int>& _BlockSizes) {
        BlockSizes = _BlockSizes;
        NumBlocks = BlockSizes.size();
        for(int size : BlockSizes) {
            if(size > 0) {
                Blocks.push_back(MatrixR::Zero(size, size));
            } else {
                Blocks.push_back(VectorR::Zero(-size));
            }
        }
    }

    inline void setIdentity() {
        for(int i = 0; i < NumBlocks; i++) {
            int size = BlockSizes[i];
            if(size > 0) {
                Blocks[i].setIdentity();
            } else {
                Blocks[i].setOnes();
            }
        }
    }

    MatrixR toDenseMatrix() const;
    DenseDiagBlock transpose() const;
    DenseDiagBlock inverse() const;
    DenseDiagBlock LLT() const;
    std::vector<Real> compute_eigenvalues() const;
    Real maxAbsCoeff() const;

    DenseDiagBlock& operator+=(const DenseDiagBlock& r);
    DenseDiagBlock& operator-=(const DenseDiagBlock& r);
    DenseDiagBlock& operator*=(const DenseDiagBlock& r);
    DenseDiagBlock& operator*=(const Real& r);

    DenseDiagBlock operator-() const;
};


inline DenseDiagBlock operator+(DenseDiagBlock lhs, const DenseDiagBlock& rhs) { return lhs += rhs; }
inline DenseDiagBlock operator-(DenseDiagBlock lhs, const DenseDiagBlock& rhs) { return lhs -= rhs; }
inline DenseDiagBlock operator*(DenseDiagBlock lhs, const DenseDiagBlock& rhs) { return lhs *= rhs; }
inline DenseDiagBlock operator*(const Real lhs, DenseDiagBlock rhs) { return rhs *= lhs; }


Real HSinner(const DenseDiagBlock& A, const DenseDiagBlock& B);





class SparseDiagBlock
{
public:
    int NumBlocks = 0;    
    std::vector<SparseMatrixR> Blocks;
    std::vector<int> BlockSizes;   // If BlockSizes[i] < 0, i-th block is diagonal matrix
    
    SparseDiagBlock() = default;
    SparseDiagBlock(const std::vector<SparseMatrixR>& _Blocks) {
        Blocks = _Blocks;
        NumBlocks = Blocks.size();
        for(auto blk : Blocks) {
            if(blk.cols() != 1) {
                assert(blk.rows() == blk.cols());
                BlockSizes.push_back( blk.rows());
            } else {
                BlockSizes.push_back(-blk.rows());
            }
        }
    }

    SparseDiagBlock(const std::vector<int>& _BlockSizes) {
        BlockSizes = _BlockSizes;
        NumBlocks = BlockSizes.size();
        for(int size : BlockSizes) {
            if(size > 0) {
                Blocks.push_back(SparseMatrixR(size, size));
            } else {
                Blocks.push_back(Eigen::SparseVector<Real>(-size));
            }
        }
    }

    MatrixR toDenseMatrix() const;

    SparseDiagBlock& operator*=(const Real& r);
    SparseDiagBlock operator-() const;

};




DenseDiagBlock& operator+=(DenseDiagBlock& lhs, const SparseDiagBlock& rhs);
DenseDiagBlock& operator-=(DenseDiagBlock& lhs, const SparseDiagBlock& rhs);
inline DenseDiagBlock operator+(DenseDiagBlock lhs, const SparseDiagBlock& rhs) {
    return lhs += rhs;
}
inline DenseDiagBlock operator-(DenseDiagBlock lhs, const SparseDiagBlock& rhs) {
    return lhs -= rhs;
}
inline SparseDiagBlock operator*(const Real& lhs, SparseDiagBlock rhs) { return rhs *= lhs; }

Real HSinner(const SparseDiagBlock& A, const DenseDiagBlock& B);


struct Options {
    unsigned int MaxIteration = 200;
    Real lambda               = 1e4;
    Real betaStar             = 0.1;
    Real betaBar              = 0.3;
    Real gamma                = 0.9;
    Real epsilon1             = 1e-30;
    Real epsilon2             = 1e-30;
    bool OUTPUT_HISTORY       = false;
};


struct Results {
    STATUS    status;
    Real  obj_primal;
    Real    obj_dual;
    VectorR        x;
    DenseDiagBlock X;
    DenseDiagBlock Y;
};


Results solve(
        const VectorR& c,
        const std::vector<SparseDiagBlock>& F,
        Options opt = Options{}
);




}