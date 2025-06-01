#include <gtest/gtest.h>

#include "core/type.hpp"
#include "math/functions.hpp"
#include "ring/Zzeta8j.hpp"

using namespace su2compiler;
using su2compiler::ring::Zzeta8j;

template <typename Real>
class Zzeta8jTest : public ::testing::Test {
protected:
    void SetUp() override {
        if constexpr(std::is_same_v<Real, mpfr::mpreal>) {
            mpfr::mpreal::set_default_prec(256);
        }
    }
};
using RealTypes = ::testing::Types<REAL_SCALAR_TYPE_LIST>;
TYPED_TEST_SUITE(Zzeta8jTest, RealTypes);


TYPED_TEST(Zzeta8jTest, Mul) {
    using Real = TypeParam;
    Integer a1 = 276;
    Integer b1 = -837;
    Integer c1 = -371;
    Integer d1 = -889;
    Integer a2 = 137;
    Integer b2 = -359;
    Integer c2 = 470;
    Integer d2 = 611;

    Integer a3 = 108;
    Integer b3 = -112;
    Integer c3 = 432;
    Integer d3 = -731;
    Integer a4 = -967;
    Integer b4 = 815;
    Integer c4 = -127;
    Integer d4 = 9;

    Zzeta8j U(a1,b1,c1,d1, a2,b2,c2,d2);
    Zzeta8j V(a3,b3,c3,d3, a3,b3,c3,d3);
    Eigen::Matrix2<Complex<Real>> U_Matrix2C = U.toMatrix2C<Real>();
    Eigen::Matrix2<Complex<Real>> V_Matrix2C = V.toMatrix2C<Real>();

    EXPECT_TRUE((U_Matrix2C * V_Matrix2C).isApprox(
            (U * V).toMatrix2C<Real>(), 
            type::epsilon<Real>()*10
    ));
}

TEST(Zzeta8jTest, Div) {
    Integer a1 = 276;
    Integer b1 = -837;
    Integer c1 = -371;
    Integer d1 = -889;
    Integer a2 = 137;
    Integer b2 = -359;
    Integer c2 = 470;
    Integer d2 = 611;

    Integer a3 = 108;
    Integer b3 = -112;
    Integer c3 = 432;
    Integer d3 = -731;
    Integer a4 = -967;
    Integer b4 = 815;
    Integer c4 = -127;
    Integer d4 = 9;

    Zzeta8j U(a1,b1,c1,d1, a2,b2,c2,d2);
    Zzeta8j V(a3,b3,c3,d3, a3,b3,c3,d3);
    Zzeta8j W = U * V;

    EXPECT_TRUE(W.leftDivisible(U));
    EXPECT_TRUE(W.rightDivisible(V));
    EXPECT_EQ(left_div(W, U), V);
    EXPECT_EQ(right_div(W, V), U);
}