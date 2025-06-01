#pragma once

#include <stdexcept>
#include <sstream>
#include <Eigen/Core>

#include "core/type.hpp"
#include "math/functions.hpp"
#include "util/time_profiler.hpp"


namespace su2compiler::math::linalg {

template <typename Real, int Rows, int Cols>
Eigen::Matrix<Real, Rows, Cols> GSO(
        const Eigen::Matrix<Real, Rows, Cols>& B
)
{
    const int M = B.rows();
    const int N = B.cols();
    Eigen::Matrix<Real, Rows, Cols> B_orth(M, N);
    auto mu = [&](int i, int j) {
        return B.col(i).dot(B_orth.col(j)) / B_orth.col(j).squaredNorm();
    }; 

    for(int i = 0; i < N; i++){
        B_orth.col(i) = B.col(i);
        for(int j = 0; j < i; j++){
            B_orth.col(i) -= mu(i, j) * B_orth.col(j);
        }
    }
    return B_orth;
}


// compute lower triangular matrix L such that A = L * L^T
template <typename Real, int Rows, int Cols>
Eigen::Matrix<Real, Rows, Cols> Cholesky(
        const Eigen::Matrix<Real, Rows, Cols>& A
)
{
    auto& prof = util::Profiler::instance();
    prof.start("Cholesky");

    if(A.rows() != A.cols()) {
        std::ostringstream oss;
        oss << "Cholesky(): "
            << "A must be square matrix, "
            << "but size of A is " << A.rows() << "x" << A.cols();
        throw std::domain_error(oss.str());
    }
    
    const int N = A.rows();
    Eigen::Matrix<Real, Rows, Cols> L(N,N);
    L.setZero();

    for(int i = 0; i < N; i++) {
        for(int j = 0; j <= i; j++) {
            Real sum(0);
            for(int k = 0; k < j; k++) {
                math::FMA(L(i,k), L(j,k), sum, sum);
            }

            if(i == j) {
                Real val = A(i,i) - sum;
                
                if(val < Real(0)) {
                    throw std::domain_error("Cholesky(): A must be symmetric positive definite matrix");
                }

                L(i,j) = math::sqrt(val);
            } else {
                L(i,j) = (A(i,j) - sum) / L(j,j);
            }
        }
    }
    
    prof.stop("Cholesky");
    return L;
}


template <typename Real, int Rows, int Cols>
Eigen::Vector<Real, Rows> SolveSystemSPD(
        const Eigen::Matrix<Real, Rows, Cols>& A,
        const Eigen::Vector<Real, Rows>& b
)
{
    auto& prof = util::Profiler::instance();
    prof.start("SolveSystemSPD");
    
    if(A.rows() != A.cols()) {
        std::ostringstream oss;
        oss << "SolveSystemSPD(): "
            << "A must be square matrix, "
            << "but size of A is " << A.rows() << "x" << A.cols();
        throw std::domain_error(oss.str());
    }

    if(A.rows() != b.rows()) {
        std::ostringstream oss;
        oss << "SolveSystemSPD(): "
            << "A and b must have same size, "
            << "but size of A is " << A.rows() << "x" << A.cols() << " "
            << "and size of b is " << b.size();
        throw std::domain_error(oss.str());
    }

    const int N = A.rows();
    Eigen::Matrix<Real, Rows, Cols> L = Cholesky(A);
    Eigen::Vector<Real, Rows> x(N);  x.setZero();
    Eigen::Vector<Real, Rows> y(N);  y.setZero();
    
    for(int i = 0; i < N; i++) {
        Real sum = b(i);
        for(int j = 0; j < i; j++) {
            sum -= L(i,j) * y(j);
        }
        y(i) = sum / L(i,i);
    }

    for(int i = N-1; i >= 0; i--) {
        Real sum = y(i);
        for(int j = i+1; j < N; j++) {
            sum -= L(j,i) * x(j);
        }
        x(i) = sum / L(i,i);
    }

    prof.stop("SolveSystemSPD");
    return x;
}


template <typename Real, int Rows, int Cols>
Eigen::MatrixX<Real> InverseSPD(
        const Eigen::Matrix<Real, Rows, Cols>& A
)
{
    auto& prof = util::Profiler::instance();
    prof.start("InberseSPD");

    if(A.rows() != A.cols()) {
        std::ostringstream oss;
        oss << "InverseSPD(): "
            << "A must be square matrix, "
            << "but size of A is " << A.rows() << "x" << A.cols();
        throw std::domain_error(oss.str());
    }

    const int N = A.rows();
    
    Eigen::Matrix<Real, Rows, Cols> L = Cholesky(A);
    
    Eigen::Matrix<Real, Rows, Cols> A_inv(N,N);
    for(int k = 0; k < N; k++) {
        Eigen::Vector<Real, Rows> b(N);  b.setZero();
        b(k) = Real(1);
        Eigen::Vector<Real, Rows> x(N);  x.setZero();
        Eigen::Vector<Real, Rows> y(N);  y.setZero();
        
        for(int i = 0; i < N; i++) {
            Real sum = b(i);
            for(int j = 0; j < i; j++) {
                sum -= L(i,j) * y(j);
            }
            y(i) = sum / L(i,i);
        }
    
        for(int i = N-1; i >= 0; i--) {
            Real sum = y(i);
            for(int j = i+1; j < N; j++) {
                sum -= L(j,i) * x(j);
            }
            x(i) = sum / L(i,i);
        }

        A_inv.col(k) = x;
    }

    prof.stop("InberseSPD");
    return A_inv; 
}


}