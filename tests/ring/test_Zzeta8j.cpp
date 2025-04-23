#include <gtest/gtest.h>

#include "type.hpp"
#include "ring/Zzeta8j.hpp"

using namespace std;
using namespace su2_compiler;
using su2_compiler::ring::Zzeta8j;


TEST(Zzeta8j, Mul) {
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
    Matrix2C U_Matrix2C = U.to_Matrix2C();
    Matrix2C V_Matrix2C = V.to_Matrix2C();

    Matrix2C D = (U * V).to_Matrix2C() - U_Matrix2C * V_Matrix2C;
    Real diff = 0.0;
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++) diff += abs(D(i,j));
    }

    const double tol = 1e-10;
    EXPECT_NEAR(static_cast<double>(diff), 0.0, tol);
}

TEST(Zzeta8j, Div) {
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