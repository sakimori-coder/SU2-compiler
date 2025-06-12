#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include "core/type.hpp"
#include "math/functions.hpp"
#include "math/linalg.hpp"
#include "util/time_profiler.hpp"

using namespace su2compiler;
using namespace su2compiler::math::linalg;



TEST(LinalgTest, GSO)
{
    Real::set_default_prec(256);

    int N = 10;
    MatrixR B = MatrixR::Random(N,N);
    
    MatrixR B_orth = GSO(B);

    Real tol = 10 * type::epsilon();

    MatrixR D = B_orth.transpose() * B_orth;
    bool isDiagonal = true;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            if(i == j) continue;
            if(math::abs(D(i,j)) < tol);
        }
    }
    EXPECT_TRUE(isDiagonal);
    
    Real val1 = B.determinant();
    Real val2 = B_orth.determinant();
    EXPECT_TRUE(math::abs(val1 - val2) < 
            std::min(math::abs(val1), math::abs(val1)) * tol);
}


TEST(LinalgTest, Cholesky)
{
    Real::set_default_prec(256);

    int N = 50;
    MatrixR Asq = MatrixR::Random(N,N);
    MatrixR A = Asq.transpose() * Asq;

    MatrixR L = Cholesky(A);

    Real tol =  10 * type::epsilon();
    EXPECT_TRUE(A.isApprox(L*L.transpose(), tol));

}


TEST(LinalgTest, SolveSystemSPD)
{
    Real::set_default_prec(256);

    int N = 50;
    MatrixR Asq = MatrixR::Random(N,N);
    MatrixR A = Asq.transpose() * Asq;
    VectorR b = VectorR::Random(N);

    VectorR x = SolveSystemSPD(A, b);

    Real tol =  1000 * type::epsilon();
    EXPECT_TRUE(b.isApprox(A*x, tol));
}


TEST(LinalgTest, InverseSPD)
{
    Real::set_default_prec(256);

    int N = 50;
    MatrixR Asq = MatrixR::Random(N,N);
    MatrixR A = Asq.transpose() * Asq;

    MatrixR A_inv = InverseSPD(A);

    Real tol =  1000 * type::epsilon();
    EXPECT_TRUE(MatrixR::Identity(N,N).isApprox(A * A_inv, tol));
}
