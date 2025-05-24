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
        if(BlockSizes[i] > 0) {
            DenseM(Eigen::seqN(START, BlockSizes[i]), Eigen::seqN(START, BlockSizes[i])) = M[i];
        } else {
            for(int j = 0; j < -BlockSizes[i]; j++) DenseM(START+j, START+j) = M[i](j);
        }
        START += std::abs(BlockSizes[i]);
    }
    return DenseM;
}

DiagBlockMatrix DiagBlockMatrix::transpose() const noexcept {
    std::vector<MatrixXR> M_transpose(NumBlocks);
    for(int i = 0; i < NumBlocks; i++) {
        if(BlockSizes[i] > 0) M_transpose[i] = M[i].transpose();
        else                  M_transpose[i] = M[i];
    }
    return DiagBlockMatrix(M_transpose);
}

DiagBlockMatrix DiagBlockMatrix::inverse() const {
    std::vector<MatrixXR> M_inv(NumBlocks);
    for(int i = 0; i < NumBlocks; i++) {
        if(BlockSizes[i] > 0) {
            M_inv[i] = M[i].inverse();
        } else {
            M_inv[i] = M[i].cwiseInverse();
        }
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
        if(BlockSizes[i] > 0){
            Eigen::LLT<MatrixXR> llt(M[i]);
            L[i] = llt.matrixL();
        } else {
            L[i] = M[i].cwiseSqrt();
        }
    }
    return DiagBlockMatrix(L);
}

std::vector<Real> DiagBlockMatrix::compute_eigenvalues() const {
    std::vector<Real> ret;
    for(int i = 0; i < NumBlocks; i++){
        if(BlockSizes[i] > 0) {
            Eigen::SelfAdjointEigenSolver<MatrixXR> es(M[i]);
            ret.insert(ret.end(), es.eigenvalues().data(), 
                                  es.eigenvalues().data() + es.eigenvalues().size());
        } else {
            ret.insert(ret.end(), M[i].data(), M[i].data() + M[i].size());
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
        if(BlockSizes[i] > 0) M[i] *= r.M[i];
        else                  M[i] = M[i].cwiseProduct(r.M[i]);
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
        if(size > 0) I.push_back(MatrixXR::Identity(size, size));
        else         I.push_back(VectorXR::Ones(-size));
    }
    return DiagBlockMatrix(I);
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
    for(int i = 0; i < A.NumBlocks; i++) {
        if(A.BlockSizes[i] > 0 && !isSymmetric(A.M[i])) return false;
    }
    return true;
}


// α = max{ λ \in [0,1] | X + λdX >= O }
Real compute_length(const DiagBlockMatrix& X, const DiagBlockMatrix& dX) {
    Real alpha;
    DiagBlockMatrix L = X.cholesky_decomposition();
    DiagBlockMatrix S = L.inverse() * dX * L.inverse().transpose();
    Real lambda_min = S.compute_eigenvalues().front();
    if(lambda_min >= 0) alpha = 1.0;
    else                alpha = min(Real(1.0), -1.0 / lambda_min);
    return alpha;
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
    int n = F[0].get_TotalSize();

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
// STEP2 : Predictor Step
//===========================================================================
    //=======================================================================
    // STEP2-1 : Search Direction
    //=======================================================================
        // Compute Xinv
        double start, end;
        start = get_time_sec();
        DiagBlockMatrix Xinv = X.inverse();
        // Compute B
        MatrixXR B(m,m);
        for(int i = 0; i < m; i++){
            DiagBlockMatrix left = Xinv * F[i+1] * Y;
            // DiagBlockMatrix left = F[i+1];
            for(int j = 0; j <= i; j++){
                DiagBlockMatrix right = F[j+1];
                B(i,j) = B(j,i) = HSinner(left, right);
            }
        }
        end = get_time_sec();
        std::cout << "1 : " << end - start << "[s]" << std::endl;
        start = get_time_sec();
        // Compute Rp
        DiagBlockMatrix Rp = -F[0] - X;
        for(int i = 0; i < m; i++) Rp += F[i+1] * x(i);
        // Compute Rc
        Real beta_p = isFeasible(epsilon1) ? 0 : betaBar;
        DiagBlockMatrix Rc_p = beta_p * mu * BlockIdentity(BlockSizes) - (X * Y);
        end = get_time_sec();
        std::cout << "2 : " << end - start << "[s]" << std::endl;
        start = get_time_sec();
        // Compute d
        VectorXR d(m);
        for(int i = 0; i < m; i++) d(i) = c(i) - HSinner(F[i+1], Y);
        // Compute r
        VectorXR r_p(m);
        DiagBlockMatrix right_p = Xinv * (Rc_p - Rp * Y);
        for(int i = 0; i < m; i++){
            DiagBlockMatrix left = F[i+1];
            r_p(i) = -d(i) + HSinner(left, right_p);
        }
        end = get_time_sec();
        std::cout << "3 : " << end - start << "[s]" << std::endl;
        start = get_time_sec();
        // Compute dx
        Eigen::LDLT<MatrixXR> ldlt(B);
        VectorXR dx_p = ldlt.solve(r_p);
        end = get_time_sec();
        std::cout << "4 : " << end - start << "[s]" << std::endl;
        // Compute dX
        DiagBlockMatrix dX_p = Rp;
        for(int i = 0; i < m; i++) dX_p += F[i+1] * dx_p(i);
        // Compute \tilde{dY}
        DiagBlockMatrix dY_tilde_p = Xinv * (Rc_p - dX_p * Y);
        // Compute dY
        DiagBlockMatrix dY_p = (dY_tilde_p + dY_tilde_p.transpose()) * (Real(1) / Real(2));        

    //=======================================================================
    // STEP2-1 : Compute Length
    //=======================================================================
        Real alpha_p_pred = compute_length(X, dX_p);
        Real alpha_d_pred = compute_length(Y, dY_p);


//===========================================================================
// STEP2 : Corrector Step
//===========================================================================
        Real beta = HSinner(X + alpha_p_pred * dX_p, Y + alpha_d_pred * dY_p) / HSinner(X, Y);
        Real beta_c;
        if(beta <= Real(1)) {
            if(isFeasible(epsilon1)) beta_c = max(betaStar, beta*beta);
            else                     beta_c = max(betaBar, beta*beta);
        } else {
                                     beta_c = Real(1);
        }
        // Compute Rc
        DiagBlockMatrix Rc = beta_c * mu * BlockIdentity(BlockSizes) - (X * Y) - (dX_p * dY_p);
        // Compute r
        VectorXR r(m);
        DiagBlockMatrix right = Xinv * (Rc - Rp * Y);
        for(int i = 0; i < m; i++){
            DiagBlockMatrix left = F[i+1];
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
        DiagBlockMatrix dY = (dY_tilde + dY_tilde.transpose()) * (Real(1) / Real(2));        

    //=======================================================================
    // STEP2-1 : Compute Length
    //=======================================================================
        Real alpha_p = compute_length(X, dX);
        Real alpha_d = compute_length(Y, dY);
        

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
            << std::setw(10)  << beta_c
            << std::endl;
        }


//===========================================================================
// STEP4 : Update
//===========================================================================
        x += gamma * alpha_p * dx;
        X += gamma * alpha_p * dX;
        Y += gamma * alpha_d * dY;
        mu = HSinner(X, Y) / Real(n);
        NumIterations++;

        // std::cout << X.compute_eigenvalues().front() << std::endl;
        // std::cout << Y.compute_eigenvalues().front() << std::endl;
    }

    return x;
}


}
}
