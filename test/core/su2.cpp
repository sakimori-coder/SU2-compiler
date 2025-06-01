#include <gtest/gtest.h>

#include <qd/qd_real.h>

#include "core/type.hpp"
#include "core/su2.hpp"

using namespace su2compiler;

template <typename Real>
class SU2Test : public ::testing::Test {
protected:
    void SetUp() override {
        if constexpr(std::is_same_v<Real, mpfr::mpreal>) {
            mpfr::mpreal::set_default_prec(256);
        }
    }
};
using RealTypes = ::testing::Types<REAL_SCALAR_TYPE_LIST>;
TYPED_TEST_SUITE(SU2Test, RealTypes);


TYPED_TEST(SU2Test, Mul) {
    using Real = TypeParam;
    SU2<Real> U = random_unitary<Real>(1234);
    SU2<Real> V = random_unitary<Real>(5678);

    Eigen::Matrix2<Complex<Real>> U_Mat = U.toEigenMatrix();
    Eigen::Matrix2<Complex<Real>> V_Mat = V.toEigenMatrix();

    EXPECT_TRUE((U_Mat * V_Mat).isApprox(
            (U*V).toEigenMatrix(),
            type::epsilon<Real>()*10
    ));
}

TYPED_TEST(SU2Test, isUnitary) {
    using Real = TypeParam;
    Real a = 100;
    Real b = 200;
    Real c = 300;
    Real d = 400;
    SU2<Real> U(a,b,c,d);
    EXPECT_FALSE(U.isUnitary());
    U.unitalize();
    EXPECT_TRUE(U.isUnitary());
}