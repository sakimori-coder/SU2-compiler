#include "probabilistic_tools/mixing_su2.hpp"

#include <Eigen/Core>

#include "type.hpp"
#include "su2.hpp"
#include "probabilistic_tools/sdp.hpp"


namespace su2compiler::mixing_su2{


MatrixXC Choi_Jamiolkowski(const MatrixXC& U) {
    int dim = U.rows();
    MatrixXC U_dag = U.adjoint();
    MatrixXC CJ_U(dim*dim, dim*dim);

    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            CJ_U(Eigen::seq(i*dim, (i+1)*dim-1), Eigen::seq(j*dim, (j+1)*dim-1)) = U.col(i) * U_dag.row(j);
        }
    }
    return CJ_U / Real(dim);
}


Matrix4R Choi_Jamiolkowski_MagicBasis(const Matrix2C& U) {
    Matrix4C CJ_U = Choi_Jamiolkowski(U);

    Complex IINV_SQRT2(Real(0), INV_SQRT2());
    Matrix4C M = Matrix4C::Zero();
    M(0,0) =  INV_SQRT2();
    M(0,3) =  INV_SQRT2();
    M(1,0) =-IINV_SQRT2;
    M(1,3) = IINV_SQRT2;
    M(2,1) =-IINV_SQRT2;
    M(2,2) =-IINV_SQRT2;
    M(3,1) =  INV_SQRT2();
    M(3,2) =- INV_SQRT2();

    Matrix4C CJMagicBasis_U = M * CJ_U * M.adjoint();
    return CJMagicBasis_U.real();
}


std::pair<Real, std::vector<Real>> optimize_distribution(const std::vector<SU2>& availableU, SU2 V) {
    const int N = availableU.size();

    std::vector<sdp::DenseMat> availableU_CJMagicBasis;
    for(auto U : availableU) availableU_CJMagicBasis.push_back(Choi_Jamiolkowski_MagicBasis(U.to_EigenMatrix()));
    sdp::DenseMat V_CJMagicBasis = Choi_Jamiolkowski_MagicBasis(V.to_EigenMatrix());

    std::vector<sdp::SparseMat> SymMat_basis(10, sdp::SparseMat(4,4));
    for(int i = 0; i < 4; i++) SymMat_basis[i].insert(i,i) = 1;
    SymMat_basis[4].insert(0,1) = SymMat_basis[4].insert(1,0) = 1;
    SymMat_basis[5].insert(0,2) = SymMat_basis[5].insert(2,0) = 1;
    SymMat_basis[6].insert(0,3) = SymMat_basis[6].insert(3,0) = 1;
    SymMat_basis[7].insert(1,2) = SymMat_basis[7].insert(2,1) = 1;
    SymMat_basis[8].insert(1,3) = SymMat_basis[8].insert(3,1) = 1;
    SymMat_basis[9].insert(2,3) = SymMat_basis[9].insert(3,2) = 1;


    // 変数の並び順は (Sの変数10個, 混合確率)
    int m = 10 + N;
    VectorXR c = VectorXR::Zero(m);
    c(0) = c(1) = c(2) = c(3) = Real(0.5);

    int NumBlocks = 4;
    std::vector<int> BlockSizes(NumBlocks);
    BlockSizes[0] = 4;      // S >= 0
    BlockSizes[1] = 4;      // S >= J(U - \sum p(x)U_x)
    BlockSizes[2] = -N;     // p(x) >= 0
    BlockSizes[3] = -1;     // \sum p(x) <= 1

    std::vector<std::vector<sdp::MatVecVar>> F_Blocks(m+1, std::vector<sdp::MatVecVar>(NumBlocks));

    // Set F0
    F_Blocks[0][0] = sdp::SparseMat(4, 4);
    F_Blocks[0][1] = V_CJMagicBasis;
    F_Blocks[0][2] = sdp::SparseVec(N);
    F_Blocks[0][3] = sdp::DenseVec(-VectorXR::Ones(1));

    // Set F1~F10
    for(int i = 1; i <= 10; i++) {
        F_Blocks[i][0] = SymMat_basis[i-1];
        F_Blocks[i][1] = SymMat_basis[i-1];
        F_Blocks[i][2] = sdp::SparseVec(N);
        F_Blocks[i][3] = sdp::DenseVec(VectorXR::Zero(1));
    }

    // Set F11~Fm
    for(int i = 0; i < N; i++) {
        F_Blocks[i+11][0] = sdp::SparseMat(4,4);
        F_Blocks[i+11][1] = availableU_CJMagicBasis[i];
        sdp::SparseVec e_i(N);
        e_i.insert(i) = 1;
        F_Blocks[i+11][2] = e_i;
        F_Blocks[i+11][3] = sdp::DenseVec(-VectorXR::Ones(1));
    }

    std::vector<sdp::DiagBlockMatrix> F;
    for(int i = 0; i <= m; i++) F.push_back(sdp::DiagBlockMatrix(F_Blocks[i]));
    Real lambda = 1e4;
    Real betaStar = 0.1;
    Real betaBar = 0.3;
    Real gamma = 0.9;
    Real epsilon1 = 1e-30;
    Real epsilon2 = 1e-30;
    int maxITERATION = 200;
    bool OUTPUT_HISTORY = true;
    VectorXR x = sdp::SDP(
        c,
        F,
        lambda,
        betaStar,
        betaBar,
        gamma,
        epsilon1,
        epsilon2,
        maxITERATION,
        OUTPUT_HISTORY
    );
    
    Real opt_distance = c.dot(x);
    std::vector<Real> opt_distribution(N);
    for(int i = 0; i < N; i++) opt_distribution[i] = x(i+10);
    
    std::cout << std::reduce(opt_distribution.begin(), opt_distribution.end()) << std::endl;

    return {opt_distance, opt_distribution};
}


}