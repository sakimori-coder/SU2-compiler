#include <gtest/gtest.h>

#include "core/type.hpp"
#include "math/functions.hpp"
#include "ring/Zzeta8.hpp"

using namespace su2compiler;
using su2compiler::ring::Zzeta8;


TEST(Zzeta8Test, Add) {
    Real::set_default_prec(256);

    Integer a1 = 8518;
    Integer b1 = 3094;
    Integer c1 = -8684;
    Integer d1 = -5115;

    Integer a2 = -3113;
    Integer b2 = 8658;
    Integer c2 = 5655;
    Integer d2 = 3751;

    Zzeta8 x(a1, b1, c1, d1);
    Zzeta8 y(a2, b2, c2, d2);
    Complex x_Complex = x.toComplex();
    Complex y_Complex = y.toComplex();

    Complex val1 = (x + y).toComplex();
    Complex val2 = x_Complex + y_Complex;

    EXPECT_TRUE(math::abs(val1 - val2) < 
            std::min(math::abs(val1), math::abs(val2)) * type::epsilon()*10);
}

TEST(Zzeta8Test, Sub) {
    Real::set_default_prec(256);

    Integer a1 = 8518;
    Integer b1 = 3094;
    Integer c1 = -8684;
    Integer d1 = -5115;

    Integer a2 = -3113;
    Integer b2 = 8658;
    Integer c2 = 5655;
    Integer d2 = 3751;

    Zzeta8 x(a1, b1, c1, d1);
    Zzeta8 y(a2, b2, c2, d2);
    Complex x_Complex = x.toComplex();
    Complex y_Complex = y.toComplex();

    Complex val1 = (x - y).toComplex();
    Complex val2 = x_Complex - y_Complex;

    EXPECT_TRUE(math::abs(val1 - val2) < 
            std::min(math::abs(val1), math::abs(val2)) * type::epsilon()*10);
}

TEST(Zzeta8Test, Mul) {
    Real::set_default_prec(256);

    Integer a1 = 8518;
    Integer b1 = 3094;
    Integer c1 = -8684;
    Integer d1 = -5115;

    Integer a2 = -3113;
    Integer b2 = 8658;
    Integer c2 = 5655;
    Integer d2 = 3751;

    a2 = 0;
    b2 = 0;
    c2 = 1;
    d2 = 0;

    Zzeta8 x(a1, b1, c1, d1);
    Zzeta8 y(a2, b2, c2, d2);
    Complex x_Complex = x.toComplex();
    Complex y_Complex = y.toComplex();

    Complex val1 = (x * y).toComplex();
    Complex val2 = x_Complex * y_Complex;

    EXPECT_TRUE(math::abs(val1 - val2) < 
            std::min(math::abs(val1), math::abs(val2)) * type::epsilon()*10);
}

TEST(Zzeta8Test, Div) {
    Integer a1 = 851;
    Integer b1 = 309;
    Integer c1 = -868;
    Integer d1 = -511;

    Integer a2 = -311;
    Integer b2 = 865;
    Integer c2 = 565;
    Integer d2 = 375;

    Zzeta8 x(a1, b1, c1, d1);
    Zzeta8 y(a2, b2, c2, d2);
    Zzeta8 z = x * y;

    EXPECT_TRUE(z.divisible(x));
    EXPECT_TRUE(z.divisible(y));
    EXPECT_EQ(z / x, y);
    EXPECT_EQ(z / y, x);
}
