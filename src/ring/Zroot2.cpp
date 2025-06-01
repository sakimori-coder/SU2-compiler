#include "ring/Zroot2.hpp"

#include <iostream>
#include <stdexcept>

#include "core/type.hpp"
#include "math/constants.hpp"

namespace su2compiler {
namespace ring {

//==============================================================================
//  implement Zroot2
//==============================================================================

//  basic properties
template <typename T>
[[nodiscard]] T Zroot2::toReal() const noexcept {
    return type::Int_to_Real<T>(a) + type::Int_to_Real<T>(b) * math::SQRT2<T>();
}
#define X(T) template [[nodiscard]] T Zroot2::toReal() const noexcept;
REAL_SCALAR_TYPES
#undef X

[[nodiscard]] bool Zroot2::divisible(const Integer& r) const{
    if(a % r != 0) return false;
    if(b % r != 0) return false;
    return true;
}

[[nodiscard]] bool Zroot2::divisible(const Zroot2& r) const{
    Zroot2 numerator = (*this) * r.conj_sqrt2();
    return numerator.divisible(r.norm_sqrt2());
}


// Zroot2 arithmetic (compound)
Zroot2& Zroot2::operator+=(const Integer& r) noexcept {
    a += r;  return *this;
}
Zroot2& Zroot2::operator+=(const Zroot2& r) noexcept {
    a += r.a;  b += r.b;  return *this;
}

Zroot2& Zroot2::operator-=(const Integer& r) noexcept {
    a -= r;  return *this;
}
Zroot2& Zroot2::operator-=(const Zroot2& r) noexcept {
    a -= r.a;  b -= r.b;  return *this;
}
Zroot2 Zroot2::operator-() const noexcept {
    return {-a, -b};
}

Zroot2& Zroot2::operator*=(const Integer& r) noexcept {
    a *= r; b *= r;  return *this;
}
Zroot2& Zroot2::operator*=(const Zroot2& r) noexcept {
    Integer na = a * r.a + 2 * b * r.b;
    Integer nb = a * r.b + b * r.a;
    a = na;  b = nb; return *this;
}

Zroot2& Zroot2::operator/=(const Integer& r) {
    if(!this->divisible(r)){
        throw std::domain_error("Zroot2::operator/= : element is not divisible by the given element");
    }
    a /= r; b /= r;  return *this;
}
Zroot2& Zroot2::operator/=(const Zroot2& r) {
    if(!this->divisible(r)){
        throw std::domain_error("Zroot2::operator/= : element is not divisible by the given element");
    }
    *this = (*this) * r.conj_sqrt2();
    Integer denominator = r.norm_sqrt2();
    *this /= denominator;
    return *this; 
}

//------------------------------------------------------------------------
//  ostream
std::ostream& operator<<(std::ostream& os, const Zroot2& x) {
    // if(x.b == 0) return os << x.a;
    // if(x.a == 0) return os << x.b << "√2";
    // return os << x.a << " + " << x.b << "√2";
    return os << "(" << x.a << ", " << x.b << ")";
}

}
}