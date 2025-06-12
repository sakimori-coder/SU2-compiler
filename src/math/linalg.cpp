#include "math/linalg.hpp"

#include <cassert>
#include <Eigen/Core>

#include "core/type.hpp"
#include "math/functions.hpp"
#include "util/time_profiler.hpp"



namespace su2compiler::math::linalg {

MatrixR GSO(const MatrixR& B)
{
    const int M = B.rows();
    const int N = B.cols();
    MatrixR B_orth(M, N);
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


MatrixR Cholesky(const MatrixR& A)
{
    auto& prof = util::Profiler::instance();
    prof.start("Cholesky");

    assert(A.rows() == A.cols());
    
    const int N = A.rows();
    MatrixR L(N,N);
    L.setZero();

    Real sum;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j <= i; j++) {
            sum = 0;
            for(int k = 0; k < j; k++) {
                math::FMA(sum, L(i,k), L(j,k), sum);
            }

            if(i == j) {
                Real val = A(i,i) - sum;

                assert(val >= Real(0));
                L(i,j) = math::sqrt(val);
            } else {
                L(i,j) = (A(i,j) - sum) / L(j,j);
            }
        }
    }
    
    prof.stop("Cholesky");
    return L;
}


VectorR SolveSystemSPD(
    const MatrixR& A,
    const VectorR& b
)
{
    auto& prof = util::Profiler::instance();
    prof.start("SolveSystemSPD");

    assert(A.rows() == A.cols());
    assert(A.rows() == b.rows());

    const int N = A.rows();
    MatrixR L = Cholesky(A);
    VectorR x(N);  x.setZero();
    VectorR y(N);  y.setZero();

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


MatrixR InverseSPD(const MatrixR& A)
{
    auto& prof = util::Profiler::instance();
    prof.start("InberseSPD");

    assert(A.rows() == A.cols());

    const int N = A.rows();
    MatrixR L = Cholesky(A);
    
    MatrixR A_inv(N,N);
    for(int k = 0; k < N; k++) {
        VectorR b(N);  b.setZero();
        b(k) = Real(1);
        VectorR x(N);  x.setZero();
        VectorR y(N);  y.setZero();
        
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