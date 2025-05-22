#include <gtest/gtest.h>

#include "type.hpp"
#include "ring/Zroot2.hpp"

using namespace std;
using namespace su2compiler;
using su2compiler::ring::Zroot2;

TEST(Zroot2, Add) {
    Integer a1 =  8695;
    Integer b1 = -3657;
    Integer a2 = -1039;
    Integer b2 =  1926;
    Zroot2 x(a1, b1);
    Zroot2 y(a2, b2);
    Real x_Real = x.to_Real();
    Real y_Real = y.to_Real();

    Real diff = (x + y).to_Real() - (x_Real + y_Real);

    const double tol = 1e-10;
    EXPECT_NEAR(static_cast<double>(diff), 0.0, tol);
}

TEST(Zroot2, Sub) {
    Integer a1 =  8695;
    Integer b1 = -3657;
    Integer a2 = -1039;
    Integer b2 =  1926;
    Zroot2 x(a1, b1);
    Zroot2 y(a2, b2);
    Real x_Real = x.to_Real();
    Real y_Real = y.to_Real();

    Real diff = (x - y).to_Real() - (x_Real - y_Real);

    const double tol = 1e-10;
    EXPECT_NEAR(static_cast<double>(diff), 0.0, tol);
}

TEST(Zroot2, Mul) {
    Integer a1 =  8695;
    Integer b1 = -3657;
    Integer a2 = -1039;
    Integer b2 =  1926;
    Zroot2 x(a1, b1);
    Zroot2 y(a2, b2);
    Real x_Real = x.to_Real();
    Real y_Real = y.to_Real();

    Real diff = (x * y).to_Real() - (x_Real * y_Real);

    const double tol = 1e-10;
    EXPECT_NEAR(static_cast<double>(diff), 0.0, tol);
}

TEST(Zroot2, Div) {
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