#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include "core/type.hpp"
#include "math/functions.hpp"
#include "math/linalg.hpp"
#include "util/time_profiler.hpp"

using namespace su2compiler;
using namespace su2compiler::math::linalg;


template <typename Real>
class LinalgTest : public ::testing::Test {
protected:
    void SetUp() override {
        if constexpr(std::is_same_v<Real, mpfr::mpreal>) {
            mpfr::mpreal::set_default_prec(256);
        }
    }
};
using RealTypes = ::testing::Types<REAL_SCALAR_TYPE_LIST>;
TYPED_TEST_SUITE(LinalgTest, RealTypes);


TYPED_TEST(LinalgTest, GSO) {
    using Real = TypeParam;

    int N = 10;
    Eigen::MatrixX<Real> B = Eigen::MatrixX<Real>::Random(N,N);
    
    Eigen::MatrixX<Real> B_orth = GSO(B);

    Real tol = 10 * type::epsilon<Real>();

    Eigen::MatrixX<Real> D = B_orth.transpose() * B_orth;
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


TYPED_TEST(LinalgTest, Cholesky) {
    using Real = TypeParam;

    int N = 50;
    Eigen::MatrixX<Real> Asq = Eigen::MatrixX<Real>::Random(N,N);
    Eigen::MatrixX<Real> A = Asq.transpose() * Asq;

    Eigen::MatrixX<Real> L = Cholesky(A);

    Real tol =  10 * type::epsilon<Real>();
    EXPECT_TRUE(A.isApprox(L*L.transpose(), tol));

}


TYPED_TEST(LinalgTest, SolveSystemSPD) {
    using Real = TypeParam;

    int N = 50;
    Eigen::MatrixX<Real> Asq = Eigen::MatrixX<Real>::Random(N,N);
    Eigen::MatrixX<Real> A = Asq.transpose() * Asq;
    Eigen::VectorX<Real> b = Eigen::VectorX<Real>::Random(N);

    Eigen::VectorX<Real> x = SolveSystemSPD(A, b);

    Real tol =  1000 * type::epsilon<Real>();
    EXPECT_TRUE(b.isApprox(A*x, tol));
}


TYPED_TEST(LinalgTest, InverseSPD) {
    using Real = TypeParam;

    int N = 50;
    Eigen::MatrixX<Real> Asq = Eigen::MatrixX<Real>::Random(N,N);
    Eigen::MatrixX<Real> A = Asq.transpose() * Asq;

    Eigen::MatrixX<Real> A_inv = InverseSPD(A);

    Real tol =  1000 * type::epsilon<Real>();
    EXPECT_TRUE(Eigen::MatrixX<Real>::Identity(N,N).isApprox(A * A_inv, tol));
}
