#include <gtest/gtest.h>

#include "core/type.hpp"
#include "core/su2.hpp"

using namespace su2compiler;


TEST(SU2Test, Mul)
{
    Real::set_default_prec(256);

    SU2 U = random_unitary(1234);
    SU2 V = random_unitary(5678);

    MatrixC U_Mat = U.toMatrixC();
    MatrixC V_Mat = V.toMatrixC();

    EXPECT_TRUE((U_Mat * V_Mat).isApprox(
                (U * V).toMatrixC(),
                type::epsilon()*10
    ));
}

TEST(SU2Test, isUnitary)
{
    Real a = 100;
    Real b = 200;
    Real c = 300;
    Real d = 400;
    SU2 U(a,b,c,d);
    EXPECT_FALSE(U.isUnitary());
    U.unitalize();
    EXPECT_TRUE(U.isUnitary());
}