#include <gtest/gtest.h>

#include "core/type.hpp"
#include "math/functions.hpp"
#include "ring/Zroot2.hpp"

using namespace su2compiler;
using su2compiler::ring::Zroot2;


TEST(Zroot2Test, Add) {
    Real::set_default_prec(256);

    Integer a1 =  8695;
    Integer b1 = -3657;
    Integer a2 = -1039;
    Integer b2 =  1926;
    Zroot2 x(a1, b1);
    Zroot2 y(a2, b2);
    Real x_Real = x.toReal();
    Real y_Real = y.toReal();

    Real val1 = (x + y).toReal();
    Real val2 = x_Real + y_Real;

    EXPECT_TRUE(math::abs(val1 - val2) < 
            std::min(math::abs(val1), math::abs(val2)) * type::epsilon()*10);
}

TEST(Zroot2Test, Sub) {
    Real::set_default_prec(256);

    Integer a1 =  8695;
    Integer b1 = -3657;
    Integer a2 = -1039;
    Integer b2 =  1926;
    Zroot2 x(a1, b1);
    Zroot2 y(a2, b2);
    Real x_Real = x.toReal();
    Real y_Real = y.toReal();

    Real val1 = (x - y).toReal();
    Real val2 = x_Real - y_Real;

    EXPECT_TRUE(math::abs(val1 - val2) < 
            std::min(math::abs(val1), math::abs(val2)) * type::epsilon()*10);
}

TEST(Zroot2Test, Mul) {
    Real::set_default_prec(256);
    
    Integer a1 =  8695;
    Integer b1 = -3657;
    Integer a2 = -1039;
    Integer b2 =  1926;
    Zroot2 x(a1, b1);
    Zroot2 y(a2, b2);
    Real x_Real = x.toReal();
    Real y_Real = y.toReal();

    Real val1 = (x * y).toReal();
    Real val2 = x_Real * y_Real;

    EXPECT_TRUE(math::abs(val1 - val2) < 
            std::min(math::abs(val1), math::abs(val2)) * type::epsilon()*10);
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