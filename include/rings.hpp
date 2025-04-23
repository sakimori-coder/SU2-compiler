#pragma once

#include <iostream>
#include <Eigen/Core>
#include "type.hpp"

namespace su2_compiler
{



/*
    Z[√2]の構造体
    a + b√2
*/
    struct ZRoot2
    {
        Integer a;
        Integer b;

        ZRoot2() : a(0), b(0) {}
        ZRoot2(Integer a, Integer b) : a(a), b(b) {}
        ZRoot2(Integer a) : a(a), b(0) {}
        ZRoot2(int a) : a(a), b(0) {}

        Integer norm() const;
        Real to_Real() const;

        bool divisibleBySqrt2() const;
        void divideBySqrt2();
    
        // 比較演算子
        bool operator<(const ZRoot2& r) const;
        bool operator<=(const ZRoot2& r) const;
        bool operator>(const ZRoot2& r) const;
        bool operator>=(const ZRoot2& r) const; 
        bool operator==(const ZRoot2& r) const;
        bool operator!=(const ZRoot2& r) const;
        // 複合代入演算子
        ZRoot2& operator+=(const ZRoot2& r);
        ZRoot2& operator+=(const Integer& r);
        ZRoot2& operator-=(const ZRoot2& r);
        ZRoot2& operator-=(const Integer& r);
        ZRoot2& operator*=(const ZRoot2& r);
        ZRoot2& operator*=(const Integer& r);
        ZRoot2& operator/=(const ZRoot2& r);
        ZRoot2& operator/=(const Integer& r);
        // 単項演算子
        ZRoot2 operator-();
    };
    
    // 算術演算
    ZRoot2 operator+(const ZRoot2& l, const ZRoot2& r);
    ZRoot2 operator+(const Integer& l, const ZRoot2& r);
    ZRoot2 operator+(const ZRoot2& l, const Integer& r);
    ZRoot2 operator-(const ZRoot2& l, const ZRoot2& r);
    ZRoot2 operator-(const Integer& l, const ZRoot2& r);
    ZRoot2 operator-(const ZRoot2& l, const Integer& r);
    ZRoot2 operator*(const ZRoot2& l, const ZRoot2& r);
    ZRoot2 operator*(const Integer& l, const ZRoot2& r);
    ZRoot2 operator*(const int& l, const ZRoot2& r);
    ZRoot2 operator*(const ZRoot2& l, const Integer& r);
    ZRoot2 operator/(const ZRoot2& l, const ZRoot2& r);
    ZRoot2 operator/(const Integer& l, const ZRoot2& r);
    ZRoot2 operator/(const ZRoot2& l, const Integer& r);    
    
    std::ostream& operator<<(std::ostream& os, const ZRoot2& x);
    
    inline ZRoot2 sqrt2_conj(const ZRoot2& x) { ZRoot2 ret(x.a, -x.b); return ret; }




/*
    Z[ω]の構造体
    a + bω + cω^2 + dω^3
*/
    struct ZOmega
    {
        Integer a;
        Integer b;
        Integer c;
        Integer d;

        ZOmega(Integer A, Integer B, Integer C, Integer D) : a(A), b(B), c(C), d(D) {}
        ZOmega(Integer A) : a(A), b(0), c(0), d(0) {}
        ZOmega() : a(0), b(0), c(0), d(0) {}

        Integer norm() const;
        Complex to_Complex() const;
        ZRoot2 real() const;
        ZRoot2 imag() const;

        bool divisibleBySqrt2() const;
        void divideBySqrt2();
    
        // 比較演算子
        bool operator<(const ZOmega& r) const;
        bool operator==(const ZOmega& r) const;
        bool operator!=(const ZOmega& r) const;
        // 複合代入演算子
        ZOmega& operator+=(const ZOmega& r);
        ZOmega& operator+=(const Integer& r);
        ZOmega& operator-=(const ZOmega& r);
        ZOmega& operator-=(const Integer& r);
        ZOmega& operator*=(const ZOmega& r);
        ZOmega& operator*=(const Integer& r);
        // 単項演算子
        ZOmega operator-() const;
    };

    // 算術演算
    ZOmega operator+(const ZOmega& l, const ZOmega& r);
    ZOmega operator+(const Integer& l, const ZOmega& r);
    ZOmega operator+(const ZOmega& l, const Integer& r);
    ZOmega operator-(const ZOmega& l, const ZOmega& r);
    ZOmega operator-(const Integer& l, const ZOmega& r);
    ZOmega operator-(const ZOmega& l, const Integer& r);
    ZOmega operator*(const ZOmega& l, const ZOmega& r);
    ZOmega operator*(const Integer& l, const ZOmega& r);
    ZOmega operator*(const ZOmega& l, const Integer& r);
    ZOmega operator/(const ZOmega& l, const ZOmega& r);
    ZOmega operator/(const Integer& l, const ZOmega& r);
    ZOmega operator/(const ZOmega& l, const Integer& r);    
    
    std::ostream& operator<<(std::ostream& os, const ZOmega& x);
    
    inline ZOmega sqrt2_conj(ZOmega x) { return {x.a, -x.b, x.c, -x.d}; }
    inline ZOmega complex_conj(ZOmega x) { return {x.a, -x.d, -x.c, -x.b}; }



    // 
    struct Quaternion_ZOmega{
        ZOmega u;
        ZOmega t;
        
        Quaternion_ZOmega() : u(1), t(0) {};
        Quaternion_ZOmega(ZOmega _u, ZOmega _t) : u(_u), t(_t) {};

        Integer norm() const;
        
        bool divisible(const Integer& r) const;
        bool divisible(const Quaternion_ZOmega& r) const;

        Quaternion_ZOmega& operator*=(const Quaternion_ZOmega& r);
        Quaternion_ZOmega& operator/=(const Quaternion_ZOmega& r);
    };

    Quaternion_ZOmega operator*(const Quaternion_ZOmega& l, const Quaternion_ZOmega& r);
    Quaternion_ZOmega operator/(const Quaternion_ZOmega& l, const Quaternion_ZOmega& r);
    Quaternion_ZOmega quaternion_conj(const Quaternion_ZOmega& x);
    Quaternion_ZOmega sqrt2_conj(const Quaternion_ZOmega& x);    
}
