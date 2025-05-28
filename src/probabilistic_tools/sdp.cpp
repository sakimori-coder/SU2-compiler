#include "probabilistic_tools/sdp.hpp"

#include <iomanip>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SparseCholesky>

#include "type.hpp"


namespace su2compiler 
{
namespace sdp
{


MatrixXR DiagBlockMatrix::to_DenseMatrix() const {
    MatrixXR TotalMat = MatrixXR::Zero(TotalSize, TotalSize);
    UINT START = 0;

    for(int i = 0; i < NumBlocks; i++){
        const UINT size = std::abs(BlockSizes[i]);

        std::visit([&](auto const& blk) {
            using T = std::decay_t<decltype(blk)>;
            if constexpr (std::is_same_v<T, DenseMat>) {
                TotalMat(Eigen::seqN(START, size), Eigen::seqN(START, size)) = blk;
            }
            else if constexpr (std::is_same_v<T, SparseMat>) {
                TotalMat(Eigen::seqN(START, size), Eigen::seqN(START, size)) = blk.toDense();
            }
            else if constexpr (std::is_same_v<T, DenseVec>) {
                TotalMat(Eigen::seqN(START, size), Eigen::seqN(START, size)) = blk.asDiagonal();
            }
            else if constexpr (std::is_same_v<T, SparseVec>) {
                TotalMat(Eigen::seqN(START, size), Eigen::seqN(START, size)) = blk.toDense().asDiagonal();
            }
        }, Blocks[i]);

        START += size;
    }
    return TotalMat;
}

DiagBlockMatrix DiagBlockMatrix::transpose() const {
    std::vector<MatVecVar> Blocks_transpose(NumBlocks);

    for(int i = 0; i < NumBlocks; i++) {
        std::visit([&](auto const& blk) {
            using T = std::decay_t<decltype(blk)>;
            if constexpr (std::is_same_v<T, DenseMat>) {
                Blocks_transpose[i] = DenseMat(blk.transpose());
            }
            else if constexpr (std::is_same_v<T, SparseMat>) {
                Blocks_transpose[i] = SparseMat(blk.transpose());
            }
            else if constexpr (std::is_same_v<T, DenseVec>) {
                Blocks_transpose[i] = blk;
            }
            else if constexpr (std::is_same_v<T, SparseVec>) {
                Blocks_transpose[i] = blk;
            }
        }, Blocks[i]);
    }
    return DiagBlockMatrix(std::move(Blocks_transpose));
}

DiagBlockMatrix DiagBlockMatrix::inverse() const {
    std::vector<MatVecVar> Blocks_inv(NumBlocks);

    for(int i = 0; i < NumBlocks; i++) {
        std::visit([&](auto const& blk) {
            using T = std::decay_t<decltype(blk)>;
            if constexpr (std::is_same_v<T, DenseMat>) {
                Blocks_inv[i] = blk.inverse().eval();
            }
            else if constexpr (std::is_same_v<T, SparseMat>) {
                Blocks_inv[i] = blk.toDense().inverse().eval();
            }
            else if constexpr (std::is_same_v<T, DenseVec>) {
                Blocks_inv[i] = blk.cwiseInverse().eval();
            }
            else if constexpr (std::is_same_v<T, SparseVec>) {
                Blocks_inv[i] = blk.toDense().cwiseInverse().eval();
            }
        }, Blocks[i]);
    }
    return DiagBlockMatrix(std::move(Blocks_inv));
}

Real DiagBlockMatrix::maxAbsCoeff() const {
    Real ret(0);

    for(int i = 0; i < NumBlocks; i++) {
        std::visit([&](auto const& blk) {
            using T = std::decay_t<decltype(blk)>;
            if constexpr (std::is_same_v<T, DenseMat>) {
                ret = max(ret, blk.cwiseAbs().maxCoeff());
            }
            else if constexpr (std::is_same_v<T, SparseMat>) {
                ret = max(ret, blk.toDense().cwiseAbs().maxCoeff());
            }
            else if constexpr (std::is_same_v<T, DenseVec>) {
                ret = max(ret, blk.cwiseAbs().maxCoeff());
            }
            else if constexpr (std::is_same_v<T, SparseVec>) {
                ret = max(ret, blk.toDense().cwiseAbs().maxCoeff());
            }
        }, Blocks[i]);
    }
    return ret;
}

DiagBlockMatrix DiagBlockMatrix::LLT() const {
    std::vector<MatVecVar> Blocks_L(NumBlocks);

    for(int i = 0; i < NumBlocks; i++) {
        std::visit([&](auto const& blk) {
            using T = std::decay_t<decltype(blk)>;
            if constexpr (std::is_same_v<T, DenseMat>) {
                Eigen::LLT<DenseMat> llt(blk);
                if(llt.info() != Eigen::Success) {
                    throw std::runtime_error("LLT() : block matrix is not SPD");
                }
                Blocks_L[i] = DenseMat(llt.matrixL());
            }
            else if constexpr (std::is_same_v<T, SparseMat>) {
                Eigen::SimplicialLLT<SparseMat> llt(blk);
                if(llt.info() != Eigen::Success) {
                    throw std::runtime_error("LLT() : block matrix is not SPD");
                }
                Blocks_L[i] = SparseMat(llt.matrixL());
            }
            else if constexpr (std::is_same_v<T, DenseVec>) {
                Blocks_L[i] = blk.cwiseSqrt().eval();
            }
            else if constexpr (std::is_same_v<T, SparseVec>) {
                Blocks_L[i] = blk.cwiseSqrt().eval();
            }
        }, Blocks[i]);
    }
    return DiagBlockMatrix(std::move(Blocks_L));
}

std::vector<Real> DiagBlockMatrix::compute_eigenvalues() const {
    std::vector<Real> ret;

    for(int i = 0; i < NumBlocks; i++){
        UINT size = std::abs(BlockSizes[i]);

        std::visit([&](auto const& blk) {
            using T = std::decay_t<decltype(blk)>;
            if constexpr (std::is_same_v<T, DenseMat>) {
                auto eigs = blk.eigenvalues().real().eval();
                ret.insert(ret.end(), eigs.data(), eigs.data() + eigs.size());
            }
            else if constexpr (std::is_same_v<T, SparseMat>) {
                auto eigs = blk.toDense().eigenvalues().real().eval();
                ret.insert(ret.end(), eigs.data(), eigs.data() + eigs.size());
            }
            else if constexpr (std::is_same_v<T, DenseVec>) {
                ret.insert(ret.end(), blk.data(), blk.data() + size);
            }
            else if constexpr (std::is_same_v<T, SparseVec>) {
                DenseVec blk_dense = blk.toDense();
                ret.insert(ret.end(), blk_dense.data(), blk_dense.data() + size);
            }
        }, Blocks[i]);
    }

    std::sort(ret.begin(), ret.end());
    return ret;
}

DiagBlockMatrix& DiagBlockMatrix::operator+=(const DiagBlockMatrix& r)
{
    if(BlockSizes != r.BlockSizes) {
        throw std::domain_error("DiagBlockMatrix::operator+= : block structure mismatch");
    }

    for(int i = 0; i < NumBlocks; i++){
        Blocks[i] = std::visit([&](const auto& blk1, const auto& blk2) -> MatVecVar {
            using T1 = std::decay_t<decltype(blk1)>;
            using T2 = std::decay_t<decltype(blk2)>;
            
            if constexpr (std::is_same_v<T1, SparseMat> && std::is_same_v<T2, SparseMat>) {
                return SparseMat(blk1 + blk2);
            }
            else if constexpr (std::is_same_v<T1, SparseVec> && std::is_same_v<T2, SparseVec>) {
                return SparseVec(blk1 + blk2);
            }
            else if constexpr (std::is_same_v<T1, SparseVec>) {
                return DenseVec(blk1.toDense() + blk2);
            }
            else if constexpr (std::is_same_v<T2, SparseVec>) {
                return DenseVec(blk1 + blk2.toDense());
            }
            else {
                if(BlockSizes[i] > 0) return DenseMat(blk1 + blk2);
                else                  return DenseVec(blk1 + blk2);
            }
        }, Blocks[i], r.Blocks[i]);
    }
    return *this;
}

DiagBlockMatrix& DiagBlockMatrix::operator-=(const DiagBlockMatrix& r)
{
    if(BlockSizes != r.BlockSizes) {
        throw std::domain_error("DiagBlockMatrix::operator-= : block structure mismatch");
    }

    for(int i = 0; i < NumBlocks; i++){
        Blocks[i] = std::visit([&](const auto& blk1, const auto& blk2) -> MatVecVar {
            using T1 = std::decay_t<decltype(blk1)>;
            using T2 = std::decay_t<decltype(blk2)>;
            
            if constexpr (std::is_same_v<T1, SparseMat> && std::is_same_v<T2, SparseMat>) {
                return SparseMat(blk1 - blk2);
            }
            else if constexpr (std::is_same_v<T1, SparseVec> && std::is_same_v<T2, SparseVec>) {
                return SparseVec(blk1 - blk2);
            }
            else if constexpr (std::is_same_v<T1, SparseVec>) {
                return DenseVec(blk1.toDense() - blk2);
            }
            else if constexpr (std::is_same_v<T2, SparseVec>) {
                return DenseVec(blk1 - blk2.toDense());
            }
            else {
                if(BlockSizes[i] > 0) return DenseMat(blk1 - blk2);
                else                  return DenseVec(blk1 - blk2);
            }
        }, Blocks[i], r.Blocks[i]);
    }
    return *this;
}

DiagBlockMatrix& DiagBlockMatrix::operator*=(const DiagBlockMatrix& r)
{
    if(BlockSizes != r.BlockSizes) {
        throw std::domain_error("DiagBlockMatrix::operator*= : operands must share the same block structure");
    }

    for(int i = 0; i < NumBlocks; i++){
        Blocks[i] = std::visit([&](const auto& blk1, const auto& blk2) -> MatVecVar {
            using T1 = std::decay_t<decltype(blk1)>;
            using T2 = std::decay_t<decltype(blk2)>;
            if constexpr (std::is_same_v<T1, SparseMat> && std::is_same_v<T2, SparseMat>) {
                return SparseMat(blk1 * blk2);
            }
            else if constexpr (std::is_same_v<T1, SparseMat>) {
                DenseMat ret_dense =  blk1 * blk2;
                return SparseMat(ret_dense.sparseView(1e-20));
            }
            else if constexpr (std::is_same_v<T2, SparseMat>) {
                DenseMat ret_dense =  blk1 * blk2;
                return SparseMat(ret_dense.sparseView(1e-20));
            }
            else if constexpr (std::is_same_v<T1, SparseVec> && std::is_same_v<T2, SparseVec>) {
                return SparseVec(blk1.cwiseProduct(blk2));
            }
            else if constexpr (std::is_same_v<T1, SparseVec>) {
                return SparseVec(blk1.cwiseProduct(blk2));
            }
            else if constexpr (std::is_same_v<T2, SparseVec>) {
                return SparseVec(blk2.cwiseProduct(blk1));
            }
            else if constexpr (std::is_same_v<T1, DenseVec>) {
                return DenseVec(blk1.cwiseProduct(blk2));
            }
            else {
                return DenseMat(blk1 * blk2);
            }
        }, Blocks[i], r.Blocks[i]);
    }

    return *this;
}

DiagBlockMatrix& DiagBlockMatrix::operator*=(const Real& r) noexcept
{
    for(int i = 0; i < NumBlocks; i++) {
        Blocks[i] = std::visit([&](const auto& blk) -> MatVecVar {
            using T = std::decay_t<decltype(blk)>;
            return T(r * blk);
        }, Blocks[i]);
    }

    return *this;
}

DiagBlockMatrix DiagBlockMatrix::operator-() const noexcept
{
    std::vector<MatVecVar> Blocks_minus(NumBlocks);

    for(int i = 0; i < NumBlocks; i++) {
        std::visit([&](const auto& blk) {
            using T = std::decay_t<decltype(blk)>;
            Blocks_minus[i] = T(-blk);
        }, Blocks[i]);
    }

    return DiagBlockMatrix(std::move(Blocks_minus));
}

std::ostream& operator<<(std::ostream& os, const DiagBlockMatrix& x)
{
    os << "DiagBlockMatrix (Num Block = "
    << x.NumBlocks << ", Total size = "
    << x.TotalSize << ")\n"
    << "-------------------------------------------\n";

    Eigen::IOFormat fmt(4, 0, "  ", "\n", " ", " ", "", "");

    for(int i = 0; i < x.NumBlocks; i++) {
        UINT size = std::abs(x.BlockSizes[i]);
        os << "[Block " << i+1 << "]" << '\n';

        std::visit([&](const auto& blk) {
            using T = std::decay_t<decltype(blk)>;
            if constexpr (std::is_same_v<T, DenseMat>) {
                os << "Type = DenseMatrix" << '\n'
                   << "Size = " << size << '\n';
                std::cout << blk.format(fmt) << std::endl;
            }
            else if constexpr (std::is_same_v<T, SparseMat>) {
                os << "Type = SparseMatrix" << '\n'
                   << "Size = " << size << '\n';
                std::cout << blk.toDense().format(fmt) << std::endl;
            }
            else if constexpr (std::is_same_v<T, DenseVec>) {
                os << "Type = DenseVector" << '\n'
                   << "Size = " << size << '\n';
                std::cout << blk.transpose().format(fmt) << std::endl;
            }
            else if constexpr (std::is_same_v<T, SparseVec>) {
                os << "Type = SparseVector" << '\n'
                   << "Size = " << size << '\n';
                std::cout << blk.toDense().transpose().format(fmt) << std::endl;
            }
        }, x.Blocks[i]);

        if(i != x.NumBlocks-1) os << std::endl;
    }
    os << "-------------------------------------------";
    return os;
}   

DiagBlockMatrix BlockIdentity(const std::vector<int>& BlockSizes)
{
    UINT NumBlocks = BlockSizes.size();
    std::vector<MatVecVar> Blocks_I(NumBlocks);
    for(int i = 0; i < NumBlocks; i++) {
        UINT size = std::abs(BlockSizes[i]);
        if(BlockSizes[i] > 0) Blocks_I[i] = DenseMat(MatrixXR::Identity(size, size));
        else                  Blocks_I[i] = DenseVec(VectorXR::Ones(size));
    }
    return DiagBlockMatrix(std::move(Blocks_I));
}

Real HSinner(const DiagBlockMatrix& A, const DiagBlockMatrix& B)
{
    if(A.BlockSizes != B.BlockSizes){
        throw std::domain_error("HSinner() : operands must share the same block structure");
    }

    Real ret(0);
    for(int i = 0; i < A.NumBlocks; i++){
        std::visit([&](const auto& blk1, const auto& blk2) {
            using T1 = std::decay_t<decltype(blk1)>;
            using T2 = std::decay_t<decltype(blk2)>;

            if constexpr (std::is_same_v<T2, SparseMat> || std::is_same_v<T2, SparseVec>) {
                ret += blk2.cwiseProduct(blk1).sum();
            }
            else {
                ret += blk1.cwiseProduct(blk2).sum();
            }
        }, A.Blocks[i], B.Blocks[i]);
    }

    return ret;
}

bool isSymmetric(const MatrixXR& A, Real tol = 1e-20) {
    return A.isApprox(A.transpose(), tol);
}

bool isSymmetric(const DiagBlockMatrix& A, Real tol = 1e-20) {
    bool ret = true;
    for(int i = 0; i < A.NumBlocks; i++) {
        std::visit([&](auto const& blk) {
            using T = std::decay_t<decltype(blk)>;
            if constexpr (std::is_same_v<T, DenseMat>) {
                if(!isSymmetric(blk)) ret = false;
            }
            else if constexpr (std::is_same_v<T, SparseMat>) {
                if(!isSymmetric(blk.toDense())) ret = false;
            }
        }, A.Blocks[i]);
    }
    return ret;
}


// α = max{ λ \in [0,1] | X + λdX >= O }
Real compute_length(const DiagBlockMatrix& X, const DiagBlockMatrix& dX) {
    Real alpha;
    DiagBlockMatrix L = X.LLT();
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
    std::vector<int> BlockSizes = F[0].get_BlockSizes();
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
        DiagBlockMatrix Xinv = X.inverse();
        // Compute B
        start = get_time_sec();
        MatrixXR B(m,m);
        for(int i = 0; i < m; i++){
            // start = get_time_sec();
            DiagBlockMatrix left = Xinv * F[i+1] * Y;
            end = get_time_sec();
            // std::cout << "HSinner : " << end - start << "[s]" << std::endl;
            // DiagBlockMatrix left = F[i+1];
            for(int j = 0; j <= i; j++){
                // DiagBlockMatrix right = F[j+1];
                B(i,j) = B(j,i) = HSinner(left, F[j+1]);
            }
        }
        // std::cout << std::setprecision(4) << B << std::endl;
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
        start = get_time_sec();
        // Compute dX
        DiagBlockMatrix dX_p = Rp;
        for(int i = 0; i < m; i++) dX_p += F[i+1] * dx_p(i);
        // Compute \tilde{dY}
        DiagBlockMatrix dY_tilde_p = Xinv * (Rc_p - dX_p * Y);
        // Compute dY
        DiagBlockMatrix dY_p = (dY_tilde_p + dY_tilde_p.transpose()) * (Real(1) / Real(2));        
        end = get_time_sec();
        std::cout << "5 : " << end - start << "[s]" << std::endl;

    //=======================================================================
    // STEP2-1 : Compute Length
    //=======================================================================
        start = get_time_sec();
        Real alpha_p_pred = compute_length(X, dX_p);
        Real alpha_d_pred = compute_length(Y, dY_p);
        end = get_time_sec();
        std::cout << "6 : " << end - start << "[s]" << std::endl;


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
