#include "rings.hpp"

#include <iostream>
#include <array>
#include <stdexcept>

#include "type.hpp"

namespace su2_compiler{



// /*
//     ZRoot2の定義
// */
//     Integer ZRoot2::norm() const { return a*a - 2*b*b; }

//     bool ZRoot2::divisibleBySqrt2() const
//     {
//         return (a % 2 == 0);
//     }

//     void ZRoot2::divideBySqrt2()
//     {
//         if (!divisibleBySqrt2()) {
//             throw std::domain_error("ZRoot2::divideBySqrt2() : not divisible by sqrt(2)");
//         }
//         Integer tmp = a;
//         a = b;
//         b = tmp / 2;
//     }

//     // 比較演算子
//     bool ZRoot2::operator<(const ZRoot2& r) const { return a < r.a || ((a == r.a) && b < r.b); }
//     bool ZRoot2::operator>(const ZRoot2& r) const { return a > r.a || ((a == r.a) && b > r.b); }
//     bool ZRoot2::operator<=(const ZRoot2& r) const { return !(*this > r); }
//     bool ZRoot2::operator>=(const ZRoot2& r) const { return !(*this < r); }
//     bool ZRoot2::operator==(const ZRoot2& r) const{ return a == r.a && b == r.b; }
//     bool ZRoot2::operator!=(const ZRoot2& r) const{ return a != r.a || b != r.b; }

//     // 複合演算子
//     ZRoot2& ZRoot2::operator+=(const ZRoot2& r)
//     {
//         a += r.a;
//         b += r.b;
//         return *this;
//     }
//     ZRoot2& ZRoot2::operator+=(const Integer& r)
//     {
//         a += r;
//         return *this;
//     }
//     ZRoot2& ZRoot2::operator-=(const ZRoot2& r)
//     {
//         a -= r.a;
//         b -= r.b;
//         return *this;
//     }
//     ZRoot2& ZRoot2::operator-=(const Integer& r)
//     {
//         a -= r;
//         return *this;
//     }
//     ZRoot2& ZRoot2::operator*=(const ZRoot2& r)
//     {
//         Integer tmp = a;
//         a = a * r.a + 2 * b * r.b;
//         b = tmp * r.b + b * r.a;
//         return *this;
//     }
//     ZRoot2& ZRoot2::operator*=(const Integer& r)
//     {
//         a *= r;
//         b *= r;
//         return *this;
//     }

//     // 単項演算子
//     ZRoot2 ZRoot2::operator-()
//     {
//         return {-a, -b};
//     }
    
//     // 算術演算子
//     ZRoot2 operator+(const ZRoot2& l, const ZRoot2& r) { return {l.a + r.a, l.b + r.b}; }
//     ZRoot2 operator+(const Integer& l, const ZRoot2& r) { return {l + r.a, r.b}; }
//     ZRoot2 operator+(const ZRoot2& l, const Integer& r) { return {l.a + r, l.b}; }
//     ZRoot2 operator-(const ZRoot2& l, const ZRoot2& r) { return {l.a - r.a, l.b - r.b}; }
//     ZRoot2 operator-(const Integer& l, const ZRoot2& r) { return {l - r.a, -r.b}; }
//     ZRoot2 operator-(const ZRoot2& l, const Integer& r) { return {l.a - r, l.b}; }
//     ZRoot2 operator*(const ZRoot2& l, const ZRoot2& r){ return {l.a*r.a + 2*l.b*r.b, l.a*r.b + l.b*r.a}; }
//     ZRoot2 operator*(const Integer& l, const ZRoot2& r) { return {l*r.a, l*r.b}; }
//     ZRoot2 operator*(const int& l, const ZRoot2& r) { return {l*r.a, l*r.b}; }
//     ZRoot2 operator*(const ZRoot2& l, const Integer& r) { return {l.a*r, l.b*r}; }

//     // 除算
//     ZRoot2& ZRoot2::operator/=(const Integer& r)
//     {
//         a /= r;
//         b /= r;
//         return *this;
//     }   
//     ZRoot2& ZRoot2::operator/=(const ZRoot2& r)
//     {
//         Integer norm = r.norm();
//         *this *= sqrt2_conj(r);
//         *this /= norm;
//         return *this;
//     }
//     ZRoot2 operator/(const ZRoot2& l, const ZRoot2& r)
//     {
//         ZRoot2 ret = l;
//         ret /= r;
//         return ret;
//     }
//     ZRoot2 operator/(const Integer& l, const ZRoot2& r)
//     {
//         ZRoot2 ret = l;
//         ret /= r;
//         return ret;
//     }
//     ZRoot2 operator/(const ZRoot2& l, const Integer& r)
//     {
//         ZRoot2 ret = l;
//         ret /= r;
//         return ret;
//     }

//     std::ostream& operator<<(std::ostream& os, const ZRoot2& x)
//     {
//         os << "(" << x.a << ", " << x.b << ")";
//         return os;    
//     }
    
    


// /*
//     ZOmegaの定義
// */
//     Integer ZOmega::norm() const
//     {
//         ZOmega x = (*this) * complex_conj(*this);
//         ZRoot2 x_ZRoot2(x.a, (x.b - x.d) / 2);
//         return x_ZRoot2.norm();
//     }

//     bool ZOmega::divisibleBySqrt2() const
//     {
//         return ((a+c) % 2 == 0 && (b+d) % 2 == 0);
//     }

//     void ZOmega::divideBySqrt2()
//     {
//         if (!divisibleBySqrt2()) {
//             throw std::domain_error("ZOmega::divideBySqrt2(): not divisible by sqrt(2)");
//         }
//         Integer _a = (b-d) / 2;
//         Integer _b = (a+c) / 2;
//         Integer _c = (b+d) / 2;
//         Integer _d = (-a+c) / 2;
//         a = _a;
//         b = _b;
//         c = _c;
//         d = _d;
//     }

//     Complex ZOmega::to_Complex() const
//     {
//         Complex ret = (Real)a
//                     + (Real)b * omega
//                     + (Real)c * omega2
//                     + (Real)d * omega3;
//         return ret;
//     }

//     ZRoot2 ZOmega::real() const
//     {
//         // Assume (b-d) is even
//         return ZRoot2(a, (b - d) / 2);
//     }

//     ZRoot2 ZOmega::imag() const
//     {
//         // Assume (b+d) is even
//         return ZRoot2(c, (b + d) / 2);
//     }

//     // 比較演算子
//     bool ZOmega::operator<(const ZOmega& r) const
//     {
//         std::array<Integer, 4> ll = {a, b, c, d};
//         std::array<Integer, 4> rr = {r.a, r.b, r.c, r.d};
//         return ll < rr;
//     }
//     bool ZOmega::operator==(const ZOmega& r) const{ return a == r.a && b == r.b && c == r.c && d == r.d; }
//     bool ZOmega::operator!=(const ZOmega& r) const{ return !(*this == r); }

//     // 複合演算子
//     ZOmega& ZOmega::operator+=(const ZOmega& r)
//     {
//         a += r.a;
//         b += r.b;
//         c += r.c;
//         d += r.d;
//         return *this;
//     }
//     ZOmega& ZOmega::operator+=(const Integer& r)
//     {
//         a += r;
//         return *this;
//     }
//     ZOmega& ZOmega::operator-=(const ZOmega& r)
//     {
//         a -= r.a;
//         b -= r.b;
//         c -= r.c;
//         d -= r.d;
//         return *this;
//     }
//     ZOmega& ZOmega::operator-=(const Integer& r)
//     {
//         a -= r;
//         return *this;
//     }
//     ZOmega& ZOmega::operator*=(const ZOmega& r)
//     {
//         Integer _a = a*r.a - b*r.d - c*r.c - d*r.b;
//         Integer _b = a*r.b + b*r.a - c*r.d - d*r.c;
//         Integer _c = a*r.c + b*r.b + c*r.a - d*r.d;
//         Integer _d = a*r.d + b*r.c + c*r.b + d*r.a;
//         this->a = _a;
//         this->b = _b;
//         this->c = _c;
//         this->d = _d;
//         return *this;
//     }
//     ZOmega& ZOmega::operator*=(const Integer& r)
//     {
//         a *= r;
//         b *= r;
//         c *= r;
//         d *= r;
//         return *this;
//     }

//     // 単項演算子
//     ZOmega ZOmega::operator-() const
//     {
//         return {-a,-b,-c,-d};
//     }

//     // 算術演算子
//     ZOmega operator+(const ZOmega& l, const ZOmega& r)
//     {
//         ZOmega ret = l;
//         ret += r;
//         return ret;
//     }
//     ZOmega operator+(const Integer& l, const ZOmega& r)
//     {
//         ZOmega ret = l;
//         ret += r;
//         return ret;
//     }
//     ZOmega operator+(const ZOmega& l, const Integer& r)
//     {
//         ZOmega ret = l;
//         ret += r;
//         return ret;
//     }
//     ZOmega operator-(const ZOmega& l, const ZOmega& r)
//     {
//         ZOmega ret = l;
//         ret -= r;
//         return ret;
//     }
//     ZOmega operator-(const Integer& l, const ZOmega& r)
//     {
//         ZOmega ret = l;
//         ret -= r;
//         return ret;
//     }
//     ZOmega operator-(const ZOmega& l, const Integer& r)
//     {
//         ZOmega ret = l;
//         ret -= r;
//         return ret;
//     }
//     ZOmega operator*(const ZOmega& l, const ZOmega& r)
//     {
//         ZOmega ret = l;
//         ret *= r;
//         return ret;
//     }
//     ZOmega operator*(const Integer& l, const ZOmega& r)
//     {
//         ZOmega ret = l;
//         ret *= r;
//         return ret;
//     }
//     ZOmega operator*(const ZOmega& l, const Integer& r)
//     {
//         ZOmega ret = l;
//         ret *= r;
//         return ret;
//     }

//     std::ostream& operator<<(std::ostream& os, const ZOmega& x)
//     {
//         os << "(" << x.a << ", " << x.b << ", " << x.c << ", " << x.d << ")";
//         return os;    
//     }



    
//     Integer Quaternion_ZOmega::norm() const
//     {
//         ZOmega x = u * complex_conj(u) + t * complex_conj(t);
//         ZRoot2 x_ZRoot2(x.a, (x.b - x.d) / 2);
//         return x_ZRoot2.norm();
//     }

//     bool Quaternion_ZOmega::divisible(const Integer& r) const
//     {
//         if(u.a % r != 0) return false;
//         if(u.b % r != 0) return false;
//         if(u.c % r != 0) return false;
//         if(u.d % r != 0) return false;
//         if(t.a % r != 0) return false;
//         if(t.b % r != 0) return false;
//         if(t.c % r != 0) return false;
//         if(t.d % r != 0) return false;
//         return true;
//     }

//     bool Quaternion_ZOmega::divisible(const Quaternion_ZOmega& r) const
//     {
//         Quaternion_ZOmega CONJ = sqrt2_conj(r * quaternion_conj(r));
//         Integer denominator = r.norm();
//         Quaternion_ZOmega numerator = (*this) * quaternion_conj(r) * CONJ;
//         return numerator.divisible(denominator);
//     }

//     Quaternion_ZOmega& Quaternion_ZOmega::operator*=(const Quaternion_ZOmega& r)
//     {
//         ZOmega _u = u * r.u - complex_conj(t) * r.t;
//         ZOmega _t = t * r.u + complex_conj(u) * r.t;
//         this->u = _u;
//         this->t = _t;
//         return *this;
//     }

//     Quaternion_ZOmega quaternion_conj(const Quaternion_ZOmega& x)
//     {
//         Quaternion_ZOmega ret(complex_conj(x.u), -x.t);
//         return ret;
//     }

//     Quaternion_ZOmega sqrt2_conj(const Quaternion_ZOmega& x)
//     {
//         Quaternion_ZOmega ret(sqrt2_conj(x.u), sqrt2_conj(x.t));
//         return ret;
//     }

//     Quaternion_ZOmega operator*(const Quaternion_ZOmega& l, const Quaternion_ZOmega& r)
//     {
//         Quaternion_ZOmega ret = l;
//         ret *= r;
//         return ret;
//     }
}