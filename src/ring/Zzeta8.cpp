#include "ring/Zzeta8.hpp"

#include <iostream>
#include <stdexcept>

#include "core/type.hpp"
#include "math/constants.hpp"
#include "ring/Zroot2.hpp"


namespace su2compiler {
namespace ring {

//==============================================================================
//  implement Zzeta8
//==============================================================================

// basic properties
template <typename T>
[[nodiscard]] std::complex<T> Zzeta8::toComplex() const {
    return type::Int_to_Real<T>(a) + type::Int_to_Real<T>(b) * math::ZETA8<T>() 
         + type::Int_to_Real<T>(c) * math::ZETA8_POW2<T>() + type::Int_to_Real<T>(d) * math::ZETA8_POW3<T>();
}
#define X(T) template [[nodiscard]] std::complex<T> Zzeta8::toComplex() const;
REAL_SCALAR_TYPES
#undef X

[[nodiscard]] bool Zzeta8::divisible(const Integer& r) const {
    if(a % r != 0) return false;
    if(b % r != 0) return false;
    if(c % r != 0) return false;
    if(d % r != 0) return false;
    return true;
}

[[nodiscard]] bool Zzeta8::divisible(const Zroot2& r) const {
    return this->divisible(Zzeta8(r));
}

[[nodiscard]] bool Zzeta8::divisible(const Zzeta8& r) const {
    Zzeta8 nume = (*this) * r.conj_complex();
    Zroot2 deno = r.norm_complex();
    nume *= deno.conj_sqrt2();
    return nume.divisible(r.norm_sqrt2());
}


// arithmetic (compound)
Zzeta8& Zzeta8::operator+=(const Integer& r) noexcept {
    a += r;  return *this;
}
Zzeta8& Zzeta8::operator+=(const Zroot2& r) noexcept {
    *this += Zzeta8(r);  return *this;
}
Zzeta8& Zzeta8::operator+=(const Zzeta8& r) noexcept {
    a += r.a;  b += r.b; c += r.c; d += r.d;  return *this;
}

Zzeta8& Zzeta8::operator-=(const Integer& r) noexcept {
    a -= r;  return *this;
}
Zzeta8& Zzeta8::operator-=(const Zroot2& r) noexcept {
    *this -= Zzeta8(r);  return *this;
}
Zzeta8& Zzeta8::operator-=(const Zzeta8& r) noexcept {
    a -= r.a;  b -= r.b; c -= r.c; d -= r.d;  return *this;
}

Zzeta8& Zzeta8::operator*=(const Integer& r) noexcept {
    a *= r; b *= r; c *= r; d *= r;  return *this;
}
Zzeta8& Zzeta8::operator*=(const Zroot2& r) noexcept {
    *this *= Zzeta8(r);  return *this;
}
Zzeta8& Zzeta8::operator*=(const Zzeta8& r) noexcept {
    Integer na = a*r.a - b*r.d - c*r.c - d*r.b;
    Integer nb = a*r.b + b*r.a - c*r.d - d*r.c;
    Integer nc = a*r.c + b*r.b + c*r.a - d*r.d;
    Integer nd = a*r.d + b*r.c + c*r.b + d*r.a;
    a = na;
    b = nb;
    c = nc;
    d = nd;
    return *this;
}

Zzeta8& Zzeta8::operator/=(const Integer& r) {
    if(!this->divisible(r)){
        throw std::domain_error("Zzeta8::operator/= : element is not divisible by the given element");
    }
    a /= r; b /= r; c /= r; d /= r;  return *this;
}
Zzeta8& Zzeta8::operator/=(const Zroot2& r) {
    if(!this->divisible(r)){
        throw std::domain_error("Zzeta8::operator/= : element is not divisible by the given element");
    }
    *this /= Zzeta8(r);  return *this;
}
Zzeta8& Zzeta8::operator/=(const Zzeta8& r) {
    if(!this->divisible(r)){
        throw std::domain_error("Zroot2::operator/= : element is not divisible by the given element");
    }
    Zzeta8 nume = (*this) * r.conj_complex();
    Zroot2 deno_conj = (r.norm_complex()).conj_sqrt2();
    nume *= Zzeta8(deno_conj.a, deno_conj.b, 0, -deno_conj.b);
    nume /= r.norm_sqrt2();
    *this = nume;
    return *this;
}

//------------------------------------------------------------------------
//  ostream
std::ostream& operator<<(std::ostream& os, const Zzeta8& x) {
    return os << "(" << x.a << ", " << x.b << ", " << x.c << ", " << x.d << ")";
}

}
}