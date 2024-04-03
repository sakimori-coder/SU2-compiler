#include "rings.hpp"
#include "type.hpp"

#include <array>

namespace SU2_Compiler{
/*
    ZRoot2の定義
*/
    ITYPE ZRoot2::norm() const { return a*a + b*b; }

    // 比較演算子
    bool ZRoot2::operator<(const ZRoot2& r) const { return a < r.a || ((a == r.a) && b < r.b); }
    bool ZRoot2::operator>(const ZRoot2& r) const { return a > r.a || ((a == r.a) && b > r.b); }
    bool ZRoot2::operator<=(const ZRoot2& r) const { return !(*this > r); }
    bool ZRoot2::operator>=(const ZRoot2& r) const { return !(*this < r); }
    bool ZRoot2::operator==(const ZRoot2& r) const{ return a == r.a && b == r.b; }
    bool ZRoot2::operator!=(const ZRoot2& r) const{ return a != r.a || b != r.b; }

    // 複合演算子
    ZRoot2& ZRoot2::operator+=(const ZRoot2& r)
    {
        a += r.a;
        b += r.b;
        return *this;
    }
    ZRoot2& ZRoot2::operator+=(const ITYPE& r)
    {
        a += r;
        return *this;
    }
    ZRoot2& ZRoot2::operator-=(const ZRoot2& r)
    {
        a -= r.a;
        b -= r.b;
        return *this;
    }
    ZRoot2& ZRoot2::operator-=(const ITYPE& r)
    {
        a -= r;
        return *this;
    }
    ZRoot2& ZRoot2::operator*=(const ZRoot2& r)
    {
        a = a * r.a + 2 * b * r.b;
        b = a * r.b + b * r.a;
        return *this;
    }
    ZRoot2& ZRoot2::operator*=(const ITYPE& r)
    {
        a *= r;
        b *= r;
        return *this;
    }

    // 単項演算子
    ZRoot2& ZRoot2::operator-()
    {
        a = -a;
        b = -b;
        return *this;
    }
    
    // 算術演算子
    ZRoot2 operator+(const ZRoot2& l, const ZRoot2& r) { return {l.a + r.a, l.b + r.b}; }
    ZRoot2 operator+(const ITYPE& l, const ZRoot2& r) { return {l + r.a, r.b}; }
    ZRoot2 operator+(const ZRoot2& l, const ITYPE& r) { return {l.a + r, l.b}; }
    ZRoot2 operator-(const ZRoot2& l, const ZRoot2& r) { return {l.a - r.a, l.b - r.b}; }
    ZRoot2 operator-(const ITYPE& l, const ZRoot2& r) { return {l - r.a, -r.b}; }
    ZRoot2 operator-(const ZRoot2& l, const ITYPE& r) { return {l.a - r, l.b}; }
    ZRoot2 operator*(const ZRoot2& l, const ZRoot2& r){ return {l.a*r.a + 2*l.b*r.b, l.a*r.b + l.b*r.a}; }
    ZRoot2 operator*(const ITYPE& l, const ZRoot2& r) { return {l*r.a, l*r.b}; }
    ZRoot2 operator*(const ZRoot2& l, const ITYPE& r) { return {l.a*r, l.b*r}; }

    // 除算
    ZRoot2& ZRoot2::operator/=(const ITYPE& r)
    {
        a /= r;
        b /= r;
        return *this;
    }   
    ZRoot2& ZRoot2::operator/=(const ZRoot2& r)
    {
        ITYPE norm = r.norm();
        *this *= r;
        *this /= norm;
        return *this;
    }
    ZRoot2 operator/(const ZRoot2& l, const ZRoot2& r)
    {
        ZRoot2 ret = l;
        ret /= r;
        return l;
    }
    ZRoot2 operator/(const ITYPE& l, const ZRoot2& r)
    {
        ZRoot2 ret = l;
        ret /= r;
        return l;
    }
    ZRoot2 operator/(const ZRoot2& l, const ITYPE& r)
    {
        ZRoot2 ret = l;
        ret /= r;
        return l;
    }

    std::ostream& operator<<(std::ostream& os, const ZRoot2& x)
    {
        os << "(" << x.a << ", " << x.b << ")";
        return os;    
    }
    
    


/*
    ZOmegaの定義
*/
    ZRoot2 ZOmega::norm() const
    {
        ZOmega ret = (*this) * conj((*this));
        return ret.real(); 
    }

    // 比較演算子
    bool ZOmega::operator<(const ZOmega& r) const
    {
        std::array<ITYPE, 4> ll = {a, b, c, d};
        std::array<ITYPE, 4> rr = {r.a, r.b, r.c, r.d};
        return ll < rr;
    }
    bool ZOmega::operator==(const ZOmega& r) const{ return a == r.a && b == r.b && c == r.c && d == r.d; }
    bool ZOmega::operator!=(const ZOmega& r) const{ return !(*this == r); }

    // 複合演算子
    ZOmega& ZOmega::operator+=(const ZOmega& r)
    {
        a += r.a;
        b += r.b;
        c += r.c;
        d += r.d;
        return *this;
    }
    ZOmega& ZOmega::operator+=(const ITYPE& r)
    {
        d += r;
        return *this;
    }
    ZOmega& ZOmega::operator-=(const ZOmega& r)
    {
        a -= r.a;
        b -= r.b;
        c -= r.c;
        d -= r.d;
        return *this;
    }
    ZOmega& ZOmega::operator-=(const ITYPE& r)
    {
        d -= r;
        return *this;
    }
    ZOmega& ZOmega::operator*=(const ZOmega& r)
    {
        ITYPE _a = a*r.d + b*r.c + c*r.b + d*r.a; 
        ITYPE _b =-a*r.a + b*r.d + c*r.c + d*r.b;
        ITYPE _c =-a*r.b - b*r.a + c*r.d + d*r.c;
        ITYPE _d =-a*r.c - b*r.b - c*r.a + d*r.d;
        this->a = _a;
        this->b = _b;
        this->c = _c;
        this->d = _d;
        return *this;
    }
    ZOmega& ZOmega::operator*=(const ITYPE& r)
    {
        a *= r;
        b *= r;
        c *= r;
        d *= r;
        return *this;
    }

    // 単項演算子
    ZOmega ZOmega::operator-() const
    {
        return {-a,-b,-c,-d};
    }

    // 算術演算子
    ZOmega operator+(const ZOmega& l, const ZOmega& r)
    {
        ZOmega ret = l;
        ret += r;
        return ret;
    }
    ZOmega operator+(const ITYPE& l, const ZOmega& r)
    {
        ZOmega ret = l;
        ret += r;
        return ret;
    }
    ZOmega operator+(const ZOmega& l, const ITYPE& r)
    {
        ZOmega ret = l;
        ret += r;
        return ret;
    }
    ZOmega operator-(const ZOmega& l, const ZOmega& r)
    {
        ZOmega ret = l;
        ret -= r;
        return ret;
    }
    ZOmega operator-(const ITYPE& l, const ZOmega& r)
    {
        ZOmega ret = l;
        ret -= r;
        return ret;
    }
    ZOmega operator-(const ZOmega& l, const ITYPE& r)
    {
        ZOmega ret = l;
        ret -= r;
        return ret;
    }
    ZOmega operator*(const ZOmega& l, const ZOmega& r)
    {
        ZOmega ret = l;
        ret *= r;
        return ret;
    }
    ZOmega operator*(const ITYPE& l, const ZOmega& r)
    {
        ZOmega ret = l;
        ret *= r;
        return ret;
    }
    ZOmega operator*(const ZOmega& l, const ITYPE& r)
    {
        ZOmega ret = l;
        ret *= r;
        return ret;
    }

    std::ostream& operator<<(std::ostream& os, const ZOmega& x)
    {
        os << "(" << x.a << ", " << x.b << ", " << x.c << ", " << x.d << ")";
        return os;    
    }


}
