#include "math/sdp.hpp"

#include <algorithm>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <core/type.hpp>
#include <math/functions.hpp>
#include <math/linalg.hpp>
#include <util/time_profiler.hpp>


namespace su2compiler::math::sdp
{

MatrixR DenseDiagBlock::toDenseMatrix() const {
    int TotalSize = 0;
    for(int size : BlockSizes) TotalSize += std::abs(size);
    MatrixR TotalMat = MatrixR::Zero(TotalSize, TotalSize);

    int START = 0;
    for(int i = 0; i < NumBlocks; i++) {
        int size = std::abs(BlockSizes[i]);
        MatrixR blk = Blocks[i];

        if(size > 0) {
            TotalMat(Eigen::seqN(START, size), Eigen::seqN(START, size)) = blk; 
        } else {
            TotalMat(Eigen::seqN(START, size), Eigen::seqN(START, size)) = blk.asDiagonal();
        }
        
        START += size;
    }
    return TotalMat;
}


DenseDiagBlock DenseDiagBlock::transpose() const {
    std::vector<MatrixR> Blocks_transpose(NumBlocks);
    
    for(int i = 0; i < NumBlocks; i++) {
        int size = BlockSizes[i];
        if(size > 0) {
            Blocks_transpose[i] = Blocks[i].transpose();
        } else {
            Blocks_transpose[i] = Blocks[i];
        }
    }
    return DenseDiagBlock(Blocks_transpose);
}


DenseDiagBlock DenseDiagBlock::inverse() const {
    std::vector<MatrixR> Blocks_inverse(NumBlocks);

    for(int i = 0; i < NumBlocks; i++) {
        int size = BlockSizes[i];
        if(size > 0) {
            Blocks_inverse[i] = Blocks[i].inverse();
        } else {
            Blocks_inverse[i] = Blocks[i].cwiseInverse();
        }
    }
    return DenseDiagBlock(Blocks_inverse);
}


DenseDiagBlock DenseDiagBlock::LLT() const {
    std::vector<MatrixR> Blocks_L(NumBlocks);

    for(int i = 0; i < NumBlocks; i++) {
        int size = BlockSizes[i];
        if(size > 0) {
            Blocks_L[i] = math::linalg::Cholesky(Blocks[i]);
        } else {
            Blocks_L[i] = Blocks[i].cwiseSqrt();
        }
    }
    return DenseDiagBlock(Blocks_L);
}


std::vector<Real> DenseDiagBlock::compute_eigenvalues() const {
    std::vector<Real> ret;

    for(int i = 0; i < NumBlocks; i++) {
        int size = BlockSizes[i];
        if(size > 0) {
            Eigen::SelfAdjointEigenSolver<MatrixR> es(Blocks[i]);
            auto eigs = es.eigenvalues();
            ret.insert(ret.end(), eigs.data(), eigs.data() + eigs.size());
        } else {
            auto eigs = Blocks[i];
            ret.insert(ret.end(), eigs.data(), eigs.data() + eigs.size());  
        }
    }
    std::sort(ret.begin(), ret.end());
    return ret;
}


Real DenseDiagBlock::maxAbsCoeff() const {
    Real ret(0);
    for(int i = 0; i < NumBlocks; i++) {
        ret = std::max(ret, Blocks[i].cwiseAbs().maxCoeff());
    }
    return ret;
}


DenseDiagBlock& DenseDiagBlock::operator+=(const DenseDiagBlock& r){
    assert(BlockSizes == r.BlockSizes);
    for(int i = 0; i < NumBlocks; i++) {
        Blocks[i] += r.Blocks[i];
    }
    return *this;
}


DenseDiagBlock& DenseDiagBlock::operator-=(const DenseDiagBlock& r){
    assert(BlockSizes == r.BlockSizes);
    for(int i = 0; i < NumBlocks; i++) {
        Blocks[i] -= r.Blocks[i];
    }
    return *this;
}


DenseDiagBlock& DenseDiagBlock::operator*=(const DenseDiagBlock& r){
    assert(BlockSizes == r.BlockSizes);
    for(int i = 0; i < NumBlocks; i++) {
        int size = BlockSizes[i];
        if(size > 0) {
            Blocks[i] *= r.Blocks[i];
        } else {
            Blocks[i] = Blocks[i].cwiseProduct(r.Blocks[i]);
        }
    }
    return *this;
}


DenseDiagBlock& DenseDiagBlock::operator*=(const Real& r){
    for(int i = 0; i < NumBlocks; i++) {
        Blocks[i] *= r;
    }
    return *this;
}


DenseDiagBlock DenseDiagBlock::operator-() const {
    std::vector<MatrixR> Blocks_minus(NumBlocks);
    for(int i = 0; i < NumBlocks; i++) {
        Blocks_minus[i] = -Blocks[i];
    }
    return DenseDiagBlock(Blocks_minus);
}


Real HSinner(const DenseDiagBlock& A, const DenseDiagBlock& B) {
    assert(A.BlockSizes == B.BlockSizes);

    Real ret(0);
    for(int i = 0; i < A.NumBlocks; i++) {
        ret += A.Blocks[i].cwiseProduct(B.Blocks[i]).sum();
    }
    return ret;
}




SparseDiagBlock& SparseDiagBlock::operator*=(const Real& r) {
    for(int i = 0; i < NumBlocks; i++) {
        Blocks[i] *= r;
    }
    return *this;
}

SparseDiagBlock SparseDiagBlock::operator-() const {
    std::vector<SparseMatrixR> Blocks_minus(NumBlocks);
    for(int i = 0; i < NumBlocks; i++) {
        Blocks_minus[i] = -Blocks[i];
    }
    return Blocks_minus;
}


DenseDiagBlock& operator+=(DenseDiagBlock& lhs, const SparseDiagBlock& rhs) {
    assert(lhs.BlockSizes == rhs.BlockSizes);

    int NumBlocks = lhs.NumBlocks;
    std::vector<int> BlockSizes = lhs.BlockSizes;
    std::vector<MatrixR> Blocks_ret(NumBlocks);

    for(int i = 0; i < NumBlocks; i++) {
        lhs.Blocks[i] += rhs.Blocks[i];
    }
    return lhs;
}


DenseDiagBlock& operator-=(DenseDiagBlock& lhs, const SparseDiagBlock& rhs) {
    assert(lhs.BlockSizes == rhs.BlockSizes);

    int NumBlocks = lhs.NumBlocks;
    std::vector<int> BlockSizes = lhs.BlockSizes;
    std::vector<MatrixR> Blocks_ret(NumBlocks);

    for(int i = 0; i < NumBlocks; i++) {
        lhs.Blocks[i] -= rhs.Blocks[i];
    }
    return lhs;
}


Real HSinner(const SparseDiagBlock& A, const DenseDiagBlock& B) {
    assert(A.BlockSizes == B.BlockSizes);

    Real ret(0);
    for(int i = 0; i < A.NumBlocks; i++) {
        ret += A.Blocks[i].cwiseProduct(B.Blocks[i]).sum();
    }
    return ret;
}


Real compute_length(const DenseDiagBlock& X, const DenseDiagBlock& dX) {
    Real alpha;
    DenseDiagBlock L = X.LLT();
    DenseDiagBlock Linv = L.inverse();
    DenseDiagBlock S = Linv * dX * Linv.transpose();
    Real lambda_min  = S.compute_eigenvalues().front();
    if(lambda_min >= 0) alpha = 1;
    else                alpha = std::min(Real(1), -Real(1) / lambda_min);
    return alpha;
}




Results solve(
    const VectorR& c,
    const std::vector<SparseDiagBlock>& F,
    Options opt)
{
    assert(!F.empty());
    assert(c.rows() + 1 == F.size());
    for(SparseDiagBlock Fi : F) {
        for(auto blk : Fi.Blocks) {
            if(blk.rows() == blk.cols()) assert(blk.isApprox(blk.transpose()));
        }
    }

    auto& prof = util::Profiler::instance();
    prof.start("sdp::solve");

    int NumBlocks = F[0].NumBlocks;
    std::vector<int> BlockSizes = F[0].BlockSizes;
    int M = c.rows();
    int N = 0;
    for(int i = 0; i < NumBlocks; i++) N += std::abs(BlockSizes[i]);
    
    DenseDiagBlock I(BlockSizes); I.setIdentity();
    DenseDiagBlock X = opt.lambda * I;
    DenseDiagBlock Y = opt.lambda * I;
    VectorR x = VectorR::Zero(M);
    Real mu = HSinner(X, Y) / N;
    int NumIterations = 0;

    // Define lambda-function for feasiblity
    auto isFeasiblePrimal = [&](Real epsilon) {
        SparseDiagBlock F0 = F[0];
        DenseDiagBlock diff_mat = X + F[0];
        for(int i = 0; i < M; i++) diff_mat -= x(i) * F[i+1];
        Real max_diff = diff_mat.maxAbsCoeff();
        if(max_diff < epsilon) return true;
        else                   return false;
    };

    auto isFeasibleDual = [&](Real epsilon) {
        Real max_diff = 0.0;
        for(int i = 0; i < M; i++) {
            max_diff = std::max(max_diff, 
                                math::abs(HSinner(F[i+1], Y) - c(i))
                               );
        }
        if(max_diff < epsilon) return true;
        else                   return false;
    };

    auto isFeasible = [&](Real epsilon) {
        if(isFeasiblePrimal(epsilon) && isFeasibleDual(epsilon)) return true;
        else                                                     return false; 
    };


    // Define lambda-function for objective value
    auto PrimalObj = [&]() {
        return c.dot(x);
    };

    auto DualObj = [&]() {
        return HSinner(F[0], Y);
    };
    
    auto RelativeGap = [&]() {
        return math::abs(PrimalObj() - DualObj()) / 
                std::max(Real(1), (math::abs(PrimalObj())) + math::abs(DualObj()) / Real(2));
    };

    
    if(opt.OUTPUT_HISTORY){
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


    while(NumIterations <= opt.MaxIteration)
    {
        if(isFeasible(opt.epsilon1) && RelativeGap() < opt.epsilon2) {
            Results ret(
                STATUS::SUCCESS,
                PrimalObj(),
                DualObj(),
                x,
                X,
                Y
            );
            prof.stop("sdp::solve");
            return ret;
        }


        DenseDiagBlock Xinv = X.inverse();
        MatrixR B = MatrixR::Zero(M,M);
        prof.start("sdp::solve compute B");
        for(int l = 0; l < NumBlocks; l++)
        {
            MatrixR Xinv_blk = Xinv.Blocks[l];
            MatrixR Y_blk    = Y.Blocks[l];
            for(int i = 0; i < M; i++)
            {
                SparseMatrixR Fi_blk = F[i+1].Blocks[l];
                MatrixR left;
                if(BlockSizes[l] > 0) left = Xinv_blk * Fi_blk * Y_blk;
                for(int j = i; j < M; j++)
                {
                    SparseMatrixR Fj_blk = F[j+1].Blocks[l];
                    Real sum(0);
                    if(BlockSizes[l] > 0) {
                        for(int outer = 0; outer < Fj_blk.outerSize(); outer++) {
                            for(SparseMatrixR::InnerIterator it(Fj_blk, outer); it; ++it) {
                                sum += left(it.row(), it.col()) * it.value();
                            }
                        }
                    } else {
                        sum = Fi_blk.cwiseProduct(Xinv_blk).cwiseProduct(Fj_blk).cwiseProduct(Y_blk).sum();
                    }

                    if(i == j) {
                        B(i,i) += sum;
                    } else {
                        B(i,j) += sum;
                        B(j,i) += sum;
                    }        
                }
            }
        }
        prof.stop("sdp::solve compute B");

        // Compute Rp
        DenseDiagBlock Rp = -X - F[0];
        for(int i = 0; i < M; i++) Rp += x(i) * F[i+1];
        // Compute Rc
        Real beta_p = isFeasible(opt.epsilon1) ? 0 : opt.betaBar;
        DenseDiagBlock Rc_p = beta_p * mu * I - (X * Y);
        // Compute d
        VectorR d(M);
        for(int i = 0; i < M; i++) d(i) = c(i) - HSinner(F[i+1], Y);
        // Compute r
        VectorR r_p(M);
        DenseDiagBlock right_p = Xinv * (Rc_p - Rp * Y);
        for(int i = 0; i < M; i++) {
            r_p(i) = -d(i) + HSinner(F[i+1], right_p);
        }
        // Compute dx_p
        VectorR dx_p = math::linalg::SolveSystemSPD(B, r_p);
        // Compute Compute dX
        DenseDiagBlock dX_p = Rp;
        for(int i = 0; i < M; i++) dX_p += dx_p(i) * F[i+1];
        // Compute \tilde{dY}
        DenseDiagBlock dY_tilde_p = Xinv * (Rc_p - dX_p * Y);
        // Compute dY
        DenseDiagBlock dY_p = Real(0.5) * (dY_tilde_p + dY_tilde_p.transpose());

        // Compute length
        Real alpha_p_pred = compute_length(X, dX_p);
        Real alpha_d_pred = compute_length(Y, dY_p);


        Real beta = HSinner(X + alpha_p_pred * dX_p, Y + alpha_d_pred * dY_p) /
                    HSinner(X, Y);
        Real beta_c;
        if(beta <= Real(1)) {
            if(isFeasible(opt.epsilon1)) beta_c = std::max(opt.betaStar, beta*beta);
            else                         beta_c = std::max(opt.betaBar,   beta*beta);
        } else {
                                         beta_c = Real(1);
        }
        // Compute Rc
        DenseDiagBlock Rc = beta_c * mu * I - (X * Y) - (dX_p * dY_p);
        // Compute r
        VectorR r(M);
        DenseDiagBlock right = Xinv * (Rc - Rp * Y);
        for(int i = 0; i < M; i++) {
            r(i) = -d(i) + HSinner(F[i+1], right);
        }
        // Compute dx
        VectorR dx = math::linalg::SolveSystemSPD(B, r);
        // Compute dX
        DenseDiagBlock dX = Rp;
        for(int i = 0; i < M; i++) dX += dx(i) * F[i+1];
        // Compute \tilde{dY}
        DenseDiagBlock dY_tilde = Xinv * (Rc - dX * Y);
        // Compute dY
        DenseDiagBlock dY = Real(0.5) * (dY_tilde + dY_tilde.transpose());

        Real alpha_p = compute_length(X, dX);
        Real alpha_d = compute_length(Y, dY);


        if(opt.OUTPUT_HISTORY){
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


        x += opt.gamma * alpha_p * dx;
        X += opt.gamma * alpha_p * dX;
        Y += opt.gamma * alpha_d * dY;
        mu = HSinner(X, Y) / N;
        NumIterations++;
    }

    Results ret(
        STATUS::FAILURE,
        PrimalObj(),
        DualObj(),
        x,
        X,
        Y
    );
    return ret;
}



}