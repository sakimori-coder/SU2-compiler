#include "ring/Zzeta8j.hpp"

#include <iostream>
#include <stdexcept>

#include "type.hpp"
#include "ring/Zroot2.hpp"
#include "ring/Zzeta8.hpp"


namespace su2compiler {
namespace ring {

//==============================================================================
//  implement Zzeta8j
//==============================================================================

// basic properties
[[nodiscard]] Matrix2C Zzeta8j::to_Matrix2C() const {
    Matrix2C U;
    U << u.to_Complex(), -t.conj_complex().to_Complex(),
         t.to_Complex(),  u.conj_complex().to_Complex();
    return U;
}

[[nodiscard]] bool Zzeta8j::divisible(const Integer& r) const{
    if(u.divisible(r) && t.divisible(r)) return true;
    return false;
}

[[nodiscard]] bool Zzeta8j::divisible(const Zroot2& r) const {
    if(u.divisible(r) && t.divisible(r)) return true;
    return false;
}

[[nodiscard]] bool Zzeta8j::divisible(const Zzeta8& r) const {
    if(u.divisible(r) && t.divisible(r)) return true;
    return false;
}

[[nodiscard]] bool Zzeta8j::leftDivisible(const Zzeta8j& r) const{
    Zzeta8j nume = r.conj_quaternion() * (*this);
    Zroot2 deno = r.norm_quaternion();
    nume = deno.conj_sqrt2() * nume;
    return nume.divisible(r.norm_sqrt2());
}

[[nodiscard]] bool Zzeta8j::rightDivisible(const Zzeta8j& r) const{
    Zzeta8j nume = (*this) * r.conj_quaternion();
    Zroot2 deno = r.norm_quaternion();
    nume = nume * deno.conj_sqrt2();
    return nume.divisible(r.norm_sqrt2());
}


// arithmetic (compound)
Zzeta8j& Zzeta8j::operator*=(const Integer& r) noexcept {
    u *= r; t *= r;  return *this;
}
Zzeta8j& Zzeta8j::operator*=(const Zroot2& r) noexcept {
    u *= r; t *= r;  return *this;
}
Zzeta8j& Zzeta8j::operator*=(const Zzeta8& r) noexcept {
    u *= r; t *= r;  return *this;
}
Zzeta8j& Zzeta8j::operator*=(const Zzeta8j& r) noexcept {
    Zzeta8 nu = u * r.u - t.conj_complex() * r.t;
    Zzeta8 nt = t * r.u + u.conj_complex() * r.t;
    u = nu;
    t = nt;
    return *this;
}

Zzeta8j& Zzeta8j::operator/=(const Integer& r) {
    if(!this->divisible(r)){
        throw std::domain_error("Zzeta8::operator/= : element is not divisible by the given element");
    }
    u /= r; t /= r;  return *this;
}
Zzeta8j& Zzeta8j::operator/=(const Zroot2& r) {
    if(!this->divisible(r)){
        throw std::domain_error("Zzeta8::operator/= : element is not divisible by the given element");
    }
    u /= r; t /= r;  return *this;
}
Zzeta8j& Zzeta8j::operator/=(const Zzeta8& r) {
    if(!this->divisible(r)){
        throw std::domain_error("Zzeta8::operator/= : element is not divisible by the given element");
    }
    u /= r; t /= r;  return *this;
}

[[nodiscard]] Zzeta8j left_div(const Zzeta8j& x, const Zzeta8j& y) {
    if(!x.leftDivisible(y)){
        throw std::domain_error("left_div : second argument does not left-divide the first exactly");
    }
    Zzeta8j nume = y.conj_quaternion() * x;
    Zroot2 deno = y.norm_quaternion();
    nume = deno.conj_sqrt2() * nume;
    return nume / y.norm_sqrt2();
}

[[nodiscard]] Zzeta8j right_div(const Zzeta8j& x, const Zzeta8j& y) {
    if(!x.rightDivisible(y)){
        throw std::domain_error("right_div : second argument does not right-divide the first exactly");
    }
    Zzeta8j nume = x * y.conj_quaternion();
    Zroot2 deno = y.norm_quaternion();
    nume = nume * deno.conj_sqrt2();
    return nume / y.norm_sqrt2();
}


//  ostream
std::ostream& operator<<(std::ostream& os, const Zzeta8j& x) {
    return os << x.u << " " << x.t;
}

}
}