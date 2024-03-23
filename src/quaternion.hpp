#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Core>

#include "type.hpp"

namespace SU2_Compiler
{
/*
    四元数クラス(a + bi + cj + dk)

    [[ a+ib, c+di],
    [-c+di, a-bi]]
    四元数a+bi+cj+dkから上のような2次元複素正方行列への写像を考えるとこれは単射環準同型となる(ハミルトン積で)。
    また、単位四元数はSU(2)に対応する。
*/
    struct quaternion
    {
        FTYPE a;
        FTYPE b;
        FTYPE c;
        FTYPE d;

        quaternion() : a(0), b(0), c(0), d(0) {};
        quaternion(FTYPE _a, FTYPE _b, FTYPE _c, FTYPE _d) : a(_a), b(_b), c(_c), d(_d) {};

        FTYPE norm() const { return sqrt(a*a + b*b + c*c + d*d); }
        FTYPE trace() const { return 2 * a; }
        Matrix2cf get_Matrix() const;
        bool is_unitary(FTYPE tol) const;
        void unitalize();

        quaternion& operator*=(const quaternion& r);
        quaternion operator*(const quaternion& r) const;
    };

    std::ostream& operator<<(std::ostream& os, const quaternion& U);
    quaternion adjoint(const quaternion& U);
    FTYPE distance(const quaternion& U, const quaternion& V);
    quaternion random_unitary(size_t seed);
    
}