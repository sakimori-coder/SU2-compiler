#include <gtest/gtest.h>

#include "core/type.hpp"
#include "math/functions.hpp"
#include "ring/Zroot2.hpp"

using namespace su2compiler;
using su2compiler::ring::Zroot2;

template <typename Real>
class Zroot2Test : public ::testing::Test {
protected:
    void SetUp() override {
        if constexpr(std::is_same_v<Real, mpfr::mpreal>) {
            mpfr::mpreal::set_default_prec(256);
        }
    }
};
using RealTypes = ::testing::Types<REAL_SCALAR_TYPE_LIST>;
TYPED_TEST_SUITE(Zroot2Test, RealTypes);


TYPED_TEST(Zroot2Test, Add) {
    using Real = TypeParam;
    Integer a1 =  8695;
    Integer b1 = -3657;
    Integer a2 = -1039;
    Integer b2 =  1926;
    Zroot2 x(a1, b1);
    Zroot2 y(a2, b2);
    Real x_Real = x.toReal<Real>();
    Real y_Real = y.toReal<Real>();

    Real val1 = (x + y).toReal<Real>();
    Real val2 = x_Real + y_Real;

    EXPECT_TRUE(math::abs(val1 - val2) < 
            std::min(math::abs(val1), math::abs(val2)) * type::epsilon<Real>()*10);
}

TYPED_TEST(Zroot2Test, Sub) {
    using Real = TypeParam;
    Integer a1 =  8695;
    Integer b1 = -3657;
    Integer a2 = -1039;
    Integer b2 =  1926;
    Zroot2 x(a1, b1);
    Zroot2 y(a2, b2);
    Real x_Real = x.toReal<Real>();
    Real y_Real = y.toReal<Real>();

    Real val1 = (x - y).toReal<Real>();
    Real val2 = x_Real - y_Real;

    EXPECT_TRUE(math::abs(val1 - val2) < 
            std::min(math::abs(val1), math::abs(val2)) * type::epsilon<Real>()*10);
}

TYPED_TEST(Zroot2Test, Mul) {
    using Real = TypeParam;
    Integer a1 =  8695;
    Integer b1 = -3657;
    Integer a2 = -1039;
    Integer b2 =  1926;
    Zroot2 x(a1, b1);
    Zroot2 y(a2, b2);
    Real x_Real = x.toReal<Real>();
    Real y_Real = y.toReal<Real>();

    Real val1 = (x * y).toReal<Real>();
    Real val2 = x_Real * y_Real;

    EXPECT_TRUE(math::abs(val1 - val2) < 
            std::min(math::abs(val1), math::abs(val2)) * type::epsilon<Real>()*10);
}

TEST(Zroot2Test, Div) {
    Integer a1 =  8695;
    Integer b1 = -3657;
    Integer a2 = -1039;
    Integer b2 =  1926;
    Zroot2 x(a1, b1);
    Zroot2 y(a2, b2);
    Zroot2 z = x * y;

    EXPECT_TRUE(z.divisible(x));
    EXPECT_TRUE(z.divisible(y));
    EXPECT_EQ(z / x, y);
    EXPECT_EQ(z / y, x);
}