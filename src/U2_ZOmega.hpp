#pragma once

#include <iostream>
#include <Eigen/Core>
#include "rings.hpp"
#include "quaternion.hpp"

namespace SU2_Compiler
{
    using Mat2ZOmega = Eigen::Matrix<ZOmega, 2, 2>;

/*
    Clifford+Tゲートで正確に記述されるU(2)の構造体
    U = [[u, -t^† * ω^l],
         [t,  u^† * ω^l]]  / √2^k
*/
    struct U2_ZOmega
    {
        ZOmega u;
        ZOmega t;
        int l;
        int k;

        U2_ZOmega() : u(1), t(0), l(0), k(0) {};   // 単位行列
        U2_ZOmega(ZOmega _u, ZOmega _t, int _l, int _k) : u(_u), t(_t), l(_l), k(_k) {};

        Mat2ZOmega get_Matrix() const;    

        U2_ZOmega& operator*=(const U2_ZOmega& r);
        U2_ZOmega operator*(const U2_ZOmega& r) const;
        U2_ZOmega operator-() const;
    };

    U2_ZOmega adjoint(const U2_ZOmega& U);
    void reduction(U2_ZOmega& U);
    quaternion to_quaternion(const U2_ZOmega& U);

    extern ZOmega omega_pow[8];
}