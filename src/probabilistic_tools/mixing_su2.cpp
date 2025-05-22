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

    std::vector<Matrix4R> availableU_CJMagicBasis;
    for(auto U : availableU) availableU_CJMagicBasis.push_back(Choi_Jamiolkowski_MagicBasis(U.to_EigenMatrix()));
    Matrix4R V_CJMagicBasis = Choi_Jamiolkowski_MagicBasis(V.to_EigenMatrix());

    std::vector<Matrix4R> SymMat_basis(10, Matrix4R::Zero());
    for(int i = 0; i < 4; i++) SymMat_basis[i](i,i) = 1;
    SymMat_basis[4](0,1) = SymMat_basis[4](1,0) = 1;
    SymMat_basis[5](0,2) = SymMat_basis[5](2,0) = 1;
    SymMat_basis[6](0,3) = SymMat_basis[6](3,0) = 1;
    SymMat_basis[7](1,2) = SymMat_basis[7](2,1) = 1;
    SymMat_basis[8](1,3) = SymMat_basis[8](3,1) = 1;
    SymMat_basis[9](2,3) = SymMat_basis[9](3,2) = 1;


    // 変数の並び順は (Sの変数10個, 混合確率)
    int m = 10 + N;
    VectorXR c = VectorXR::Zero(m);
    c(0) = c(1) = c(2) = c(3) = Real(0.5);

    int NumBlocks = 3 + N;
    std::vector<int> BlockSizes(NumBlocks);
    BlockSizes[0] = 4;                                // S >= 0
    BlockSizes[1] = 4;                                // S >= J(U - \sum p(x)U_x)
    for(int i = 0; i < N; i++) BlockSizes[2+i] = 1;   // p(x) >= 0
    BlockSizes[2+N] = 1;                              // \sum p(x) <= 1

    std::vector<std::vector<MatrixXR>> Fmat(m+1, std::vector<MatrixXR>(NumBlocks));

    // Set F0
    Fmat[0][0] = Matrix4R::Zero();
    Fmat[0][1] = V_CJMagicBasis;
    for(int i = 0; i < N; i++) Fmat[0][2+i] = Matrix1R::Zero();
    Fmat[0][2+N] = -Matrix1R::Identity();

    // Set F1~F10
    for(int i = 1; i <= 10; i++) {
        Fmat[i][0] = SymMat_basis[i-1];
        Fmat[i][1] = SymMat_basis[i-1];
        for(int j = 0; j < N; j++) Fmat[i][2+j] = Matrix1R::Zero();
        Fmat[i][2+N] = Matrix1R::Zero();
    }

    // Set F11~Fm
    for(int i = 0; i < N; i++) {
        Fmat[i+11][0] = Matrix4R::Zero();
        Fmat[i+11][1] = availableU_CJMagicBasis[i];
        for(int j = 0; j < N; j++) Fmat[i+11][2+j] = Matrix1R::Zero();
        Fmat[i+11][2+i] = Matrix1R::Identity();
        Fmat[i+11][2+N] = -Matrix1R::Identity();
    }

    std::vector<sdp::DiagBlockMatrix> F;
    for(int i = 0; i <= m; i++) F.push_back(sdp::DiagBlockMatrix(Fmat[i]));
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