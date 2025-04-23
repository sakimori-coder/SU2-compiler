#include <gtest/gtest.h>

#include "type.hpp"
#include "ring/Zzeta8.hpp"

using namespace std;
using namespace su2_compiler;
using su2_compiler::ring::Zzeta8;


TEST(Zzeta8, Add) {
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
    Complex x_Complex = x.to_Complex();
    Complex y_Complex = y.to_Complex();

    Real diff = abs((x + y).to_Complex() - (x_Complex + y_Complex));

    const double tol = 1e-10;
    EXPECT_NEAR(static_cast<double>(diff), 0.0, tol);
}

TEST(Zzeta8, Sub) {
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
    Complex x_Complex = x.to_Complex();
    Complex y_Complex = y.to_Complex();

    Real diff = abs((x - y).to_Complex() - (x_Complex - y_Complex));

    const double tol = 1e-10;
    EXPECT_NEAR(static_cast<double>(diff), 0.0, tol);
}

TEST(Zzeta8, Mul) {
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
    Complex x_Complex = x.to_Complex();
    Complex y_Complex = y.to_Complex();

    Real diff = abs((x * y).to_Complex() - (x_Complex * y_Complex));

    const double tol = 1e-10;
    EXPECT_NEAR(static_cast<double>(diff), 0.0, tol);
}

TEST(Zzeta8, Div) {
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
