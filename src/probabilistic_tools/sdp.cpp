#include "probabilistic_tools/sdp.hpp"

#include <iomanip>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include "type.hpp"


namespace su2compiler 
{
namespace sdp
{


MatrixXR DiagBlockMatrix::to_DenseMatrix() const noexcept {
    MatrixXR DenseM = MatrixXR::Zero(TotalSize, TotalSize);
    int START = 0;
    for(int i = 0; i < NumBlocks; i++){
        DenseM(Eigen::seqN(START, BlockSizes[i]), Eigen::seqN(START, BlockSizes[i])) = M[i];
        START += BlockSizes[i];
    }
    return DenseM;
}

DiagBlockMatrix DiagBlockMatrix::transpose() const noexcept {
    std::vector<MatrixXR> M_transpose(NumBlocks);
    for(int i = 0; i < NumBlocks; i++) {
        M_transpose[i] = M[i].transpose();
    }
    return DiagBlockMatrix(M_transpose);
}

DiagBlockMatrix DiagBlockMatrix::inverse() const {
    std::vector<MatrixXR> M_inv(NumBlocks);
    for(int i = 0; i < NumBlocks; i++) {
        if(M[i].cols() == 1) {
            M_inv[i] = Matrix1R::Constant(Real(1.0) / M[i](0,0));
        }
        else M_inv[i] = M[i].inverse();
    }
    return DiagBlockMatrix(M_inv);
}

Real DiagBlockMatrix::maxAbsCoeff() const noexcept {
    Real ret = 0.0;
    for(int i = 0; i < NumBlocks; i++) {
        ret = max(ret, M[i].array().abs().maxCoeff());
    }
    return ret;
}

DiagBlockMatrix DiagBlockMatrix::cholesky_decomposition() const {
    std::vector<MatrixXR> L(NumBlocks);
    for(int i = 0; i < NumBlocks; i++){
        if(M[i].cols() == 1){
            L[i] = Matrix1R::Constant(M[i](0,0));
        }else {
            Eigen::LLT<MatrixXR> llt(M[i]);
            L[i] = llt.matrixL();
        }
    }
    return DiagBlockMatrix(L);
}

std::vector<Real> DiagBlockMatrix::compute_eigenvalues() const {
    std::vector<Real> ret;
    for(int i = 0; i < NumBlocks; i++){
        if(M[i].cols() == 1) {
            ret.push_back(Real(1.0));
        } else {
            Eigen::SelfAdjointEigenSolver<MatrixXR> es(M[i]);
            ret.insert(ret.end(), es.eigenvalues().data(), 
                                  es.eigenvalues().data() + es.eigenvalues().size());
        }
    }
    std::sort(ret.begin(), ret.end());
    return ret;
}

DiagBlockMatrix& DiagBlockMatrix::operator+=(const DiagBlockMatrix& r)
{
    if(BlockSizes != r.BlockSizes) {
        throw std::domain_error("DiagBlockMatrix::operator+= : operands must share the same block structure");
    }

    for(int i = 0; i < NumBlocks; i++){
        M[i] += r.M[i];
    }
    return *this;
}

DiagBlockMatrix& DiagBlockMatrix::operator-=(const DiagBlockMatrix& r)
{
    if(BlockSizes != r.BlockSizes) {
        throw std::domain_error("DiagBlockMatrix::operator-= : operands must share the same block structure");
    }

    for(int i = 0; i < NumBlocks; i++){
        M[i] -= r.M[i];
    }
    return *this;
}

DiagBlockMatrix& DiagBlockMatrix::operator*=(const DiagBlockMatrix& r)
{
    if(BlockSizes != r.BlockSizes) {
        throw std::domain_error("DiagBlockMatrix::operator*= : operands must share the same block structure");
    }

    for(int i = 0; i < NumBlocks; i++){
        M[i] *= r.M[i];
    }
    return *this;
}

DiagBlockMatrix& DiagBlockMatrix::operator*=(const Real& r) noexcept
{
    for(int i = 0; i < NumBlocks; i++){
        M[i] *= r;
    }
    return *this;
}

DiagBlockMatrix DiagBlockMatrix::operator-() const noexcept
{
    std::vector<MatrixXR> minusM = M;
    for(int i = 0; i < NumBlocks; i++) minusM[i] = -minusM[i];
    return DiagBlockMatrix(minusM);
}

std::ostream& operator<<(std::ostream& os, const DiagBlockMatrix& x)
{
    os << "DiagBlockMatrix (Num Block = "
    << x.NumBlocks << ", Total size = "
    << x.TotalSize << ")\n"
    << "-------------------------------------------\n";

    Eigen::IOFormat fmt(4, 0, "  ", "\n", " ", " ", "", "");

    for(int i = 0; i < x.NumBlocks; i++) {
        const auto& mat = x.M[i];
        os << "[Block " << i+1 << "] size = "
           << mat.rows() << " x " << mat.cols() << '\n'
           << mat.format(fmt) << std::endl;
        if(i != x.NumBlocks-1) os << std::endl;
    }
    os << "-------------------------------------------";
    return os;
    return os;
}   

DiagBlockMatrix BlockIdentity(const std::vector<int>& structure)
{
    std::vector<MatrixXR> I;
    for(int size : structure){
        I.push_back(MatrixXR::Identity(size, size));
    }
    return I;
}

Real HSinner(const DiagBlockMatrix& A, const DiagBlockMatrix& B)
{
    if(A.BlockSizes != B.BlockSizes){
        throw std::domain_error("HSinner() : operands must share the same block structure");
    }
    Real ret = 0.0;
    for(int i = 0; i < A.NumBlocks; i++){
        const MatrixXR mat1 = A.M[i];
        const MatrixXR mat2 = B.M[i];
        ret += mat1.cwiseProduct(mat2).sum();
    }
    return ret;
}

bool isSymmetric(const MatrixXR& A, Real tol = 1e-20) {
    return A.isApprox(A.transpose(), tol);
}

bool isSymmetric(const DiagBlockMatrix& A, Real tol = 1e-20) {
    for(int i = 0; i < A.NumBlocks; i++) if(!isSymmetric(A.M[i])) return false;
    return true;
}


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
    bool OUTPUT_HISTORY = false) 
{
    int m = c.rows();   // F = {F_0, F_1, ... , F_m}
    std::vector<int> BlockSizes = F[0].structure();
    int n = std::reduce(BlockSizes.begin(), BlockSizes.end());

    for(int i = 0; i <= m; i++) {
        if(!isSymmetric(F[i])) {
            throw std::domain_error("SDP() : F must be a symmetric matrix");
        }
    }


//===========================================================================
// STEP0 : Initialization
//===========================================================================
    DiagBlockMatrix X = lambda * BlockIdentity(BlockSizes);
    DiagBlockMatrix Y = lambda * BlockIdentity(BlockSizes);
    VectorXR x = VectorXR::Zero(m);
    Real mu = HSinner(X, Y) / n;
    int NumIterations = 0;

    
    // Define Lambda-function for feasiblity
    auto isFeasiblePrimal = [&](Real epsilon) {
        DiagBlockMatrix Diff_mat = X + F[0];
        for(int i = 0; i < m; i++) Diff_mat -= F[i+1] * x(i);
        Real diff = Diff_mat.maxAbsCoeff();
        if(diff < epsilon) return true;
        else               return false;
    };

    auto isFeasibleDual = [&](Real epsilon) {
        Real diff = 0.0;
        for(int i = 0; i < m; i++) {
            diff = max(diff, abs(HSinner(F[i+1], Y) - c(i)));
        }
        if(diff < epsilon) return true;
        else               return false; 
    };

    auto isFeasible = [&](Real epsilon) {
        if(isFeasiblePrimal(epsilon) && isFeasibleDual(epsilon)) return true;
        else                                                     return false;
    };


    // Define Lambda-function for objective value
    auto PrimalObj = [&]() {
        return c.dot(x);
    };

    auto DualObj = [&]() {
        return HSinner(F[0], Y);
    };

    auto RelativeGap = [&]() {
        return abs(PrimalObj() - DualObj()) / 
                    max(Real(1.0), (abs(PrimalObj()) + abs(DualObj())) / 2.0);
    };


    if(OUTPUT_HISTORY){
        std::cout << "-------------------------------------------------------------------------" << std::endl;
        std::cout << "                             SDP parameters                              " << std::endl; 
        std::cout << "-------------------------------------------------------------------------" << std::endl;
        std::cout << std::setw(8)  << "      mu"
        << std::setw(12) << "objP"
        << std::setw(10) << "objD"
        << std::setw(10) << "rgap"
        << std::setw(12) << "alphaP"
        << std::setw(10) << "alphaD"
        << std::setw(8) << "beta"
        << std::endl;
    }

    while(NumIterations <= maxITERATION)
    {
//===========================================================================
// STEP1 : Cheking Convergence
//===========================================================================
        if(isFeasible(epsilon1) && RelativeGap() < epsilon2) break;



//===========================================================================
// STEP2 : Search Direction
//===========================================================================
        // Compute Xinv
        DiagBlockMatrix Xinv = X.inverse();
        // Compute B
        MatrixXR B(m,m);
        for(int i = 0; i < m; i++){
            for(int j = 0; j < m; j++){
                DiagBlockMatrix left = Xinv * F[i+1] * Y;
                DiagBlockMatrix right = F[j+1];
                B(i,j) = HSinner(left, right);
            }
        }
        // Compute Rp
        DiagBlockMatrix Rp = -F[0] - X;
        for(int i = 0; i < m; i++) Rp += F[i+1] * x(i);
        // Compute Rc
        Real beta = isFeasible(epsilon1) ? betaStar : betaBar;
        DiagBlockMatrix Rc = beta * mu * BlockIdentity(BlockSizes) - (X * Y);
        // Compute d
        VectorXR d(m);
        for(int i = 0; i < m; i++) d(i) = c(i) - HSinner(F[i+1], Y);
        // Compute r
        VectorXR r(m);
        for(int i = 0; i < m; i++){
            DiagBlockMatrix left = F[i+1];
            DiagBlockMatrix right = Xinv * (Rc - Rp * Y);
            r(i) = -d(i) + HSinner(left, right);
        }
        // Compute dx
        VectorXR dx = B.inverse() * r;
        // Compute dX
        DiagBlockMatrix dX = Rp;
        for(int i = 0; i < m; i++) dX += F[i+1] * dx(i);
        // Compute \tilde{dY}
        DiagBlockMatrix dY_tilde = Xinv * (Rc - dX * Y);
        // Compute dY
        DiagBlockMatrix dY = (dY_tilde + dY_tilde.transpose()) * Real(0.5);
        

//===========================================================================
// STEP3 : Step Length
//===========================================================================
        // Compute α_p
        Real alpha_p;
        DiagBlockMatrix L_X = X.cholesky_decomposition();
        DiagBlockMatrix S_X = L_X.inverse() * dX * L_X.inverse().transpose();
        Real lambda_min_X = S_X.compute_eigenvalues().front();
        if(lambda_min_X >= 0) alpha_p = 1.0;
        else                  alpha_p = min(Real(1.0), -1.0 / lambda_min_X);
        // Compute α_d
        Real alpha_d;
        DiagBlockMatrix L_Y = Y.cholesky_decomposition();
        DiagBlockMatrix S_Y = L_Y.inverse() * dY * L_Y.inverse().transpose();
        Real lambda_min_Y = S_Y.compute_eigenvalues().front();
        if(lambda_min_Y >= 0) alpha_d = 1.0;
        else                  alpha_d = min(Real(1.0), -1.0 / lambda_min_Y);


//===========================================================================
// STEP3' : Output history
//===========================================================================
        if(OUTPUT_HISTORY){
            std::cout << std::setw(3)  << NumIterations
            << "  " << std::scientific << std::setprecision(1)
            << std::setw(8)  << mu
            << std::setw(10) << PrimalObj()
            << std::setw(10) << DualObj()
            << std::setw(10) << RelativeGap()
            << std::setw(10) << alpha_p
            << std::setw(10) << alpha_d
            << std::setw(10)  << beta
            << std::endl;
        }


//===========================================================================
// STEP4 : Update
//===========================================================================
        x += gamma * alpha_p * dx;
        X += gamma * alpha_p * dX;
        Y += gamma * alpha_d * dY;
        mu = HSinner(X, Y) / n;
        NumIterations++;
    }

    return x;
}


}
}
