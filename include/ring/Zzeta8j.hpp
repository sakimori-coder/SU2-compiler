#pragma once

#include <iostream>

#include "type.hpp"
#include "Zroot2.hpp"
#include "Zzeta8.hpp"

namespace su2_compiler {
namespace ring {

//==============================================================================
//  struct Zzeta8j
//==============================================================================
//  Lightweight value type for an algebraic integer in
//  the mixed cyclotomic–quaternionic ring
//
//      ℤ[ζ₈ , j ζ₈]  ≔  { a₁ + b₁ζ₈ + c₁ζ₈² + d₁ζ₈³ + j(a₂ + b₂ζ₈ + c₂ζ₈² + d₂ζ₈³) | aₖ, bₖ, cₖ, dₖ ∈ ℤ } 
//                       { u + jt  |  u,t ∈ ℤ[ζ₈] }
//
//  where  ζ₈ = e^{iπ/4}  is a primitive 8‑th root of unity and  j  is the
//  quaternion unit satisfying  j² = −1  and  ij = −ji.
//
//  Storage
//  -------
//      Zzeta8 u   coefficient of 1
//      Zzeta8 t   coefficient of j
//
//  Example
//  -------
//      using su2_compiler::ring::Zzeta8j;
//
//      Zzeta8j u(1, 1, 0, 0,            // 1 + ζ₈
//                0, 0, 0, 0);           // (no j‑part)
//
//      Zzeta8j v(0, 0, 0, 0,            // 0
//                0, 2, 0, -1);          // 2 j ζ₈¹ − 1 j ζ₈³
//
//      Zzeta8j w = u * v;               // ring multiplication
//------------------------------------------------------------------------------
struct Zzeta8j
{
    Zzeta8 u;
    Zzeta8 t;
    
    //------------------------------------------------------------------------
    //  constructor
    constexpr Zzeta8j() noexcept = default;
    constexpr Zzeta8j(Integer a1, Integer b1, Integer c1, Integer d1,
                      Integer a2, Integer b2, Integer c2, Integer d2) noexcept : u({a1,b1,c1,d1}), t({a2,b2,c2,d2}) {}
    constexpr Zzeta8j(Zzeta8 _u, Zzeta8 _t) noexcept : u(_u), t(_t) {}
    constexpr explicit Zzeta8j(Zroot2 r) noexcept : u(r), t(0) {}
    constexpr explicit Zzeta8j(Zzeta8 _u) noexcept : u(_u), t(0) {}


    //------------------------------------------------------------------------
    //  basic properties
    [[nodiscard]] Zroot2 norm_quaternion() const noexcept {
        return u.norm_complex() + t.norm_complex();
    }

    [[nodiscard]] Integer norm_sqrt2() const noexcept {
        return this->norm_quaternion().norm_sqrt2();
    }

    [[nodiscard]] constexpr Zzeta8j conj_quaternion() const noexcept {
        return {u.conj_complex(), -t};
    }

    [[nodiscard]] constexpr Zzeta8j conj_sqrt2() const noexcept {
        return {u.conj_sqrt2(), t.conj_sqrt2()};
    }

    [[nodiscard]] Matrix2C to_Matrix2C() const;

    [[nodiscard]] bool divisible(const Integer& r) const;
    [[nodiscard]] bool divisible(const Zroot2& r) const;
    [[nodiscard]] bool divisible(const Zzeta8& r) const;
    [[nodiscard]] bool leftDivisible(const Zzeta8j& r) const;
    [[nodiscard]] bool rightDivisible(const Zzeta8j& r) const;



    //------------------------------------------------------------------------
    //  arithmetic (compound)
    constexpr Zzeta8j operator-() const noexcept { return {-u, -t}; }

    Zzeta8j& operator*=(const Integer& r) noexcept;
    Zzeta8j& operator*=(const Zroot2& r) noexcept;
    Zzeta8j& operator*=(const Zzeta8& r) noexcept;
    Zzeta8j& operator*=(const Zzeta8j& r) noexcept;

    Zzeta8j& operator/=(const Integer& r);
    Zzeta8j& operator/=(const Zroot2& r);
    Zzeta8j& operator/=(const Zzeta8& r);

    auto operator<=>(const Zzeta8j&) const = default;
};

//------------------------------------------------------------------------
//  arithmetic (binary)
[[nodiscard]] inline Zzeta8j operator*(Zzeta8j lhs, const Zzeta8j& rhs) noexcept { return lhs *= rhs; }
[[nodiscard]] inline Zzeta8j operator*(Zzeta8j lhs, const Integer& rhs) noexcept { return lhs *= rhs; }
[[nodiscard]] inline Zzeta8j operator*(const Integer& lhs, Zzeta8j rhs) noexcept { return rhs *= lhs; }
[[nodiscard]] inline Zzeta8j operator*(Zzeta8j lhs, const Zroot2& rhs) noexcept { return lhs *= rhs; }
[[nodiscard]] inline Zzeta8j operator*(const Zroot2& lhs, Zzeta8j rhs) noexcept { return rhs *= lhs; }
[[nodiscard]] inline Zzeta8j operator*(Zzeta8j lhs, const Zzeta8& rhs) noexcept { return lhs *= rhs; }
[[nodiscard]] inline Zzeta8j operator*(const Zzeta8 lhs, Zzeta8j rhs) noexcept { return rhs *= lhs; }

[[nodiscard]] inline Zzeta8j operator/(Zzeta8j lhs, const Integer& rhs) { return lhs /= rhs; }
[[nodiscard]] inline Zzeta8j operator/(Zzeta8j lhs, const Zroot2& rhs) noexcept { return lhs /= rhs; }
[[nodiscard]] inline Zzeta8j operator/(Zzeta8j lhs, const Zzeta8& rhs) noexcept { return lhs /= rhs; }
// compute  y⁻¹x
[[nodiscard]] Zzeta8j left_div(const Zzeta8j& x, const Zzeta8j& y);
// compute xy⁻¹
[[nodiscard]] Zzeta8j right_div(const Zzeta8j& x, const Zzeta8j& y);

//------------------------------------------------------------------------
//  ostream
std::ostream& operator<<(std::ostream& os, const Zzeta8j& x);

}
}