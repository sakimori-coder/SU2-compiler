#include "core/mix_su2.hpp"

#include <cassert>
#include <vector>
#include <Eigen/Core>

#include "core/type.hpp"
#include "core/su2.hpp"
#include "math/constants.hpp"


namespace su2compiler
{


MatrixC Choi_Jamiolkowski(const MatrixC& U) {
    assert(U.rows() == U.cols());

    int dim = U.rows();
    MatrixC U_dag = U.adjoint();
    MatrixC CJ_U(dim*dim, dim*dim);

    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            CJ_U(Eigen::seq(i*dim, (i+1)*dim-1),
                 Eigen::seq(j*dim, (j+1)*dim-1)) = U.col(i) * U_dag.row(j);
        }
    }
    return CJ_U / Real(dim);
}


MatrixR Choi_Jamiolkowski_MagicBasis(const MatrixC& U) {
    assert(U.rows() == 2 && U.cols() == 2);
    
    MatrixC CJ_U = Choi_Jamiolkowski(U);

    using math::INV_SQRT2;
    Complex IINV_SQRT2(Real(0), INV_SQRT2());
    MatrixC M = MatrixC::Zero(4,4);
    M(0,0) =  INV_SQRT2();
    M(0,3) =  INV_SQRT2();
    M(1,0) =-IINV_SQRT2;
    M(1,3) = IINV_SQRT2;
    M(2,1) =-IINV_SQRT2;
    M(2,2) =-IINV_SQRT2;
    M(3,1) =  INV_SQRT2();
    M(3,2) =- INV_SQRT2();

    MatrixC CJMagicBasis_U = M * CJ_U * M.adjoint();
    return CJMagicBasis_U.real();
}



using namespace math::sdp;
Results MixSU2::compute_optimal_prob(
    const SU2& targetV,
    Options opt)
{
    Real::set_default_prec(256);

    const int N = availableU.size();
    
    std::vector<SparseMatrixR> availableU_CJMB;
    for(auto U : availableU) {
        availableU_CJMB.push_back(Choi_Jamiolkowski_MagicBasis(U.toMatrixC()).sparseView());
    }
    SparseMatrixR V_CJMB = Choi_Jamiolkowski_MagicBasis(targetV.toMatrixC()).sparseView();

    std::vector<SparseMatrixR> SymMat_basis(10, SparseMatrixR(4,4));
    for(int i = 0; i < 4; i++) SymMat_basis[i].insert(i,i) = 1;
    SymMat_basis[4].insert(0,1) = SymMat_basis[4].insert(1,0) = 1;
    SymMat_basis[5].insert(0,2) = SymMat_basis[5].insert(2,0) = 1;
    SymMat_basis[6].insert(0,3) = SymMat_basis[6].insert(3,0) = 1;
    SymMat_basis[7].insert(1,2) = SymMat_basis[7].insert(2,1) = 1;
    SymMat_basis[8].insert(1,3) = SymMat_basis[8].insert(3,1) = 1;
    SymMat_basis[9].insert(2,3) = SymMat_basis[9].insert(3,2) = 1;

    int M = 10 + N;
    VectorR c = VectorR::Zero(M);
    c(0) = c(1) = c(2) = c(3) = Real(0.5);

    int NumBlocks = 4;
    std::vector<int> BlockSizes(NumBlocks);
    BlockSizes[0] = 4;      // S >= 0
    BlockSizes[1] = 4;      // S >= J(U - \sum p(x)U_x)
    BlockSizes[2] = -N;     // p(x) >= 0
    BlockSizes[3] = -1;     // \sum p(x) <= 1

    std::vector<SparseDiagBlock> F(M+1);

    SparseMatrixR ScalarOne(1,1);  ScalarOne.insert(0,0) = 1;
    // Set F0
    std::vector<SparseMatrixR> F0_Blocks(NumBlocks);
    F0_Blocks[0] = SparseMatrixR(4,4);
    F0_Blocks[1] = V_CJMB;
    F0_Blocks[2] = SparseMatrixR(N, 1);
    F0_Blocks[3] = -ScalarOne;
    F[0] = SparseDiagBlock(F0_Blocks);

    // Set F1~F10
    for(int i = 1; i <= 10; i++) {
        std::vector<SparseMatrixR> Fi_Blocks(NumBlocks);
        Fi_Blocks[0] = SymMat_basis[i-1];
        Fi_Blocks[1] = SymMat_basis[i-1];
        Fi_Blocks[2] = SparseMatrixR(N, 1);
        Fi_Blocks[3] = SparseMatrixR(1,1);
        F[i] = SparseDiagBlock(Fi_Blocks);
    }

    // Set F11~F_M
    for(int i = 0; i < N; i++) {
        std::vector<SparseMatrixR> Fi_Blocks(NumBlocks);
        Fi_Blocks[0] = SparseMatrixR(4,4);
        Fi_Blocks[1] = availableU_CJMB[i];
        SparseMatrixR e_i(N,1); e_i.insert(i,0) = 1;
        Fi_Blocks[2] = e_i;
        Fi_Blocks[3] = -ScalarOne;
        F[i+11] = SparseDiagBlock(Fi_Blocks);
    }

    Results ret = solve(c, F, opt);
    
    
    
    return ret;
}



}