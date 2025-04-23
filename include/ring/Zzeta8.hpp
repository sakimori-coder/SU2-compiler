#pragma once

#include <iostream>

#include "type.hpp"
#include "Zroot2.hpp"


namespace su2_compiler {
namespace ring {

//==============================================================================
//  struct Zzeta8
//==============================================================================
//  Lightweight value type for an algebraic integer in the cyclotomic ring
//  ℤ[ζ₈] ≔ { a + bζ₈ + cζ₈² + dζ₈³ | a,b,c,d ∈ ℤ } , where  ζ₈ = e^{iπ/4}  
//  is a primitive 8‑th root of unity.
//
//  Storage
//  -------
//      Integer a   coefficient of 1
//      Integer b   coefficient of ζ₈
//      Integer c   coefficient of ζ₈²
//      Integer d   coefficient of ζ₈³
//
//  Example
//  -------
//      using su2_compiler::ring::Zzeta8;
//
//      Zzeta8 u(1, 0, 1, 0);            // 1 + ζ₈²
//      Zzeta8 v(0, 1, 0, -1);           // ζ₈ – ζ₈³
//      Zzeta8 w = u * v;                // ring multiplication
//
//      w = (1 + ζ₈²)(ζ₈ – ζ₈³)
//        = ζ₈ − ζ₈³ + ζ₈³ − ζ₈⁵
//        = ζ₈ + ζ₈  (since ζ₈⁴ = −1)
//        = 2ζ₈
//------------------------------------------------------------------------------
struct Zzeta8
{
    Integer a = 0;
    Integer b = 0;
    Integer c = 0;
    Integer d = 0;
    

    //------------------------------------------------------------------------
    //  constructor
    constexpr Zzeta8() noexcept = default;
    constexpr explicit Zzeta8(Integer _a) noexcept : a(_a), b(0), c(0), d(0) {}
    constexpr explicit Zzeta8(Zroot2 r) noexcept : a(r.a), b(r.b), c(0), d(-r.b) {}
    constexpr Zzeta8(Integer _a, Integer _b, Integer _c, Integer _d) noexcept : a(_a), b(_b), c(_c), d(_d) {}
    
    //------------------------------------------------------------------------
    //  basic properties
    [[nodiscard]] constexpr Zroot2 norm_complex() const noexcept {
        return {a*a + b*b + c*c + d*d, a*b - a*d + c*b + c*d};
    }

    [[nodiscard]] Integer norm_sqrt2() const noexcept {
        return this->norm_complex().norm_sqrt2();
    }

    [[nodiscard]] constexpr Zzeta8 conj_complex() const noexcept {
        return {a, -d, -c, -b};
    }

    [[nodiscard]] constexpr Zzeta8 conj_sqrt2() const noexcept {
        return {a, -b, c, -d};
    }

    [[nodiscard]] Complex to_Complex() const;

    [[nodiscard]] bool divisible(const Integer& r) const;
    [[nodiscard]] bool divisible(const Zroot2& r) const;
    [[nodiscard]] bool divisible(const Zzeta8& r) const;



    //------------------------------------------------------------------------
    //  arithmetic (compound)
    Zzeta8& operator+=(const Integer& r) noexcept;
    Zzeta8& operator+=(const Zroot2& r) noexcept;
    Zzeta8& operator+=(const Zzeta8& r) noexcept;

    Zzeta8& operator-=(const Integer& r) noexcept;
    Zzeta8& operator-=(const Zroot2& r) noexcept;
    Zzeta8& operator-=(const Zzeta8& r) noexcept;

    Zzeta8& operator*=(const Integer& r) noexcept;
    Zzeta8& operator*=(const Zroot2& r) noexcept;
    Zzeta8& operator*=(const Zzeta8& r) noexcept;

    Zzeta8& operator/=(const Integer& r);
    Zzeta8& operator/=(const Zroot2& r);
    Zzeta8& operator/=(const Zzeta8& r);

    constexpr Zzeta8 operator-() const noexcept { return {-a,-b,-c,-d}; }

    auto operator<=>(const Zzeta8&) const = default;
};

//------------------------------------------------------------------------
//  arithmetic (binary)
[[nodiscard]] inline Zzeta8 operator+(Zzeta8 lhs, const Zzeta8& rhs) noexcept { return lhs += rhs; }
[[nodiscard]] inline Zzeta8 operator+(Zzeta8 lhs, const Integer& rhs) noexcept { return lhs += rhs; }
[[nodiscard]] inline Zzeta8 operator+(const Integer& lhs, Zzeta8 rhs) noexcept { return rhs += lhs; }
[[nodiscard]] inline Zzeta8 operator+(Zzeta8 lhs, const Zroot2& rhs) noexcept { return lhs += rhs; }
[[nodiscard]] inline Zzeta8 operator+(const Zroot2& lhs, Zzeta8 rhs) noexcept { return rhs += lhs; }

[[nodiscard]] inline Zzeta8 operator-(Zzeta8 lhs, const Zzeta8& rhs) noexcept { return lhs -= rhs; }
[[nodiscard]] inline Zzeta8 operator-(Zzeta8 lhs, const Integer& rhs) noexcept { return lhs -= rhs; }
[[nodiscard]] inline Zzeta8 operator-(const Integer& lhs, Zzeta8 rhs) noexcept { return rhs -= lhs; }
[[nodiscard]] inline Zzeta8 operator-(Zzeta8 lhs, const Zroot2& rhs) noexcept { return lhs -= rhs; }
[[nodiscard]] inline Zzeta8 operator-(const Zroot2& lhs, Zzeta8 rhs) noexcept { return rhs -= lhs; }

[[nodiscard]] inline Zzeta8 operator*(Zzeta8 lhs, const Zzeta8& rhs) noexcept { return lhs *= rhs; }
[[nodiscard]] inline Zzeta8 operator*(Zzeta8 lhs, const Integer& rhs) noexcept { return lhs *= rhs; }
[[nodiscard]] inline Zzeta8 operator*(const Integer& lhs, Zzeta8 rhs) noexcept { return rhs *= lhs; }
[[nodiscard]] inline Zzeta8 operator*(Zzeta8 lhs, const Zroot2& rhs) noexcept { return lhs *= rhs; }
[[nodiscard]] inline Zzeta8 operator*(const Zroot2& lhs, Zzeta8 rhs) noexcept { return rhs *= lhs; }

[[nodiscard]] inline Zzeta8 operator/(Zzeta8 lhs, const Zzeta8& rhs) { return lhs /= rhs; }
[[nodiscard]] inline Zzeta8 operator/(Zzeta8 lhs, const Integer& rhs) { return lhs /= rhs; }
[[nodiscard]] inline Zzeta8 operator/(Zzeta8 lhs, const Zroot2& rhs) { return lhs /= rhs; }

//------------------------------------------------------------------------
//  ostream
std::ostream& operator<<(std::ostream& os, const Zzeta8& x);

}
}