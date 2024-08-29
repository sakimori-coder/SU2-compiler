#pragma once

#include <iostream>
#include <tbb/concurrent_hash_map.h>
#include "type.hpp"

namespace SU2_Compiler
{
/*
    Z[√2]の構造体
    a + b√2
*/
    struct ZRoot2
    {
        ITYPE a;
        ITYPE b;

        ZRoot2() : a(0), b(0) {}
        ZRoot2(ITYPE a, ITYPE b) : a(a), b(b) {}
        ZRoot2(ITYPE a) : a(a), b(0) {}

        ITYPE norm() const;
    
        // 比較演算子
        bool operator<(const ZRoot2& r) const;
        bool operator<=(const ZRoot2& r) const;
        bool operator>(const ZRoot2& r) const;
        bool operator>=(const ZRoot2& r) const; 
        bool operator==(const ZRoot2& r) const;
        bool operator!=(const ZRoot2& r) const;
        // 複合代入演算子
        ZRoot2& operator+=(const ZRoot2& r);
        ZRoot2& operator+=(const ITYPE& r);
        ZRoot2& operator-=(const ZRoot2& r);
        ZRoot2& operator-=(const ITYPE& r);
        ZRoot2& operator*=(const ZRoot2& r);
        ZRoot2& operator*=(const ITYPE& r);
        ZRoot2& operator/=(const ZRoot2& r);
        ZRoot2& operator/=(const ITYPE& r);
        // 単項演算子
        ZRoot2 operator-();
    };
    
    // 算術演算
    ZRoot2 operator+(const ZRoot2& l, const ZRoot2& r);
    ZRoot2 operator+(const ITYPE& l, const ZRoot2& r);
    ZRoot2 operator+(const ZRoot2& l, const ITYPE& r);
    ZRoot2 operator-(const ZRoot2& l, const ZRoot2& r);
    ZRoot2 operator-(const ITYPE& l, const ZRoot2& r);
    ZRoot2 operator-(const ZRoot2& l, const ITYPE& r);
    ZRoot2 operator*(const ZRoot2& l, const ZRoot2& r);
    ZRoot2 operator*(const ITYPE& l, const ZRoot2& r);
    ZRoot2 operator*(const ZRoot2& l, const ITYPE& r);
    ZRoot2 operator/(const ZRoot2& l, const ZRoot2& r);
    ZRoot2 operator/(const ITYPE& l, const ZRoot2& r);
    ZRoot2 operator/(const ZRoot2& l, const ITYPE& r);    
    
    std::ostream& operator<<(std::ostream& os, const ZRoot2& x);
    
    inline FTYPE ZRoot2_to_FTYPE(ZRoot2& x) { return (FTYPE)x.a + sqrt2*(FTYPE)x.b; }
    inline ZRoot2 conj(const ZRoot2& x) { ZRoot2 ret(x.a, -x.b); return ret; }




/*
    Z[ω]の構造体
    aω^3 + bω^2 + cω + d
*/
    struct ZOmega
    {
        ITYPE a;
        ITYPE b;
        ITYPE c;
        ITYPE d;

        ZOmega(ITYPE A, ITYPE B, ITYPE C, ITYPE D) : a(A), b(B), c(C), d(D) {}
        ZOmega(ITYPE D) : a(0), b(0), c(0), d(D) {}
        ZOmega() : a(0), b(0), c(0), d(0) {}

        ZRoot2 norm() const;
        // aω^3+bω^2+cω+d = α + iβ (α,β \in Z[√2]を仮定 <=> c+aとc-aが2で割り切れる)
        ZRoot2 real() const { return {d, (c-a) / 2}; }   // αを返す
        ZRoot2 imag() const { return {b, (c+a) / 2}; }   // βを返す
    
        // 比較演算子
        bool operator<(const ZOmega& r) const;
        bool operator==(const ZOmega& r) const;
        bool operator!=(const ZOmega& r) const;
        // 複合代入演算子
        ZOmega& operator+=(const ZOmega& r);
        ZOmega& operator+=(const ITYPE& r);
        ZOmega& operator-=(const ZOmega& r);
        ZOmega& operator-=(const ITYPE& r);
        ZOmega& operator*=(const ZOmega& r);
        ZOmega& operator*=(const ITYPE& r);
        ZOmega& operator/=(const ZOmega& r);
        ZOmega& operator/=(const ITYPE& r);
        // 単項演算子
        ZOmega operator-() const;
    };

    // 算術演算
    ZOmega operator+(const ZOmega& l, const ZOmega& r);
    ZOmega operator+(const ITYPE& l, const ZOmega& r);
    ZOmega operator+(const ZOmega& l, const ITYPE& r);
    ZOmega operator-(const ZOmega& l, const ZOmega& r);
    ZOmega operator-(const ITYPE& l, const ZOmega& r);
    ZOmega operator-(const ZOmega& l, const ITYPE& r);
    ZOmega operator*(const ZOmega& l, const ZOmega& r);
    ZOmega operator*(const ITYPE& l, const ZOmega& r);
    ZOmega operator*(const ZOmega& l, const ITYPE& r);
    ZOmega operator/(const ZOmega& l, const ZOmega& r);
    ZOmega operator/(const ITYPE& l, const ZOmega& r);
    ZOmega operator/(const ZOmega& l, const ITYPE& r);    
    
    std::ostream& operator<<(std::ostream& os, const ZOmega& x);
    
    inline CTYPE ZOmega_to_FTYPE(const ZOmega& x) { return (FTYPE)x.a*omega3 + (FTYPE)x.b*omega2 + (FTYPE)x.c*omega + (FTYPE)x.d; }
    inline ZOmega conj(ZOmega x) { return {-x.c, -x.b, -x.a, x.d}; }

    struct ZRoot2_hash
    {
        
        std::size_t operator()(const ZRoot2& key) const {
            // ハッシュ関数を定義（簡単な例）
            std::size_t seed = std::hash<ITYPE>{}(key.a);
            return std::hash<ITYPE>{}(key.b) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }

        std::size_t hash(const ZRoot2& key) const {
            // ハッシュ関数を定義（簡単な例）
            std::size_t seed = std::hash<ITYPE>{}(key.a);
            return std::hash<ITYPE>{}(key.b) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }

        bool equal(const ZRoot2& x, const ZRoot2& y) const{
            return x == y;
        }
    };
}



namespace std {
    template <>
    struct hash<SU2_Compiler::ZRoot2> {
        private:
            const std::hash<SU2_Compiler::ITYPE> ITYPE_hash; 

        public:
            hash() : ITYPE_hash() {}
            std::size_t operator()(const SU2_Compiler::ZRoot2& key) const {
                // ハッシュ関数を定義（簡単な例）
                std::size_t seed = ITYPE_hash(key.a);
                return ITYPE_hash(key.b) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            }
    };
}
