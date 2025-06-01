#pragma once

#include <iostream>

#include "core/type.hpp"

namespace su2compiler {
namespace ring {

//==============================================================================
//  struct Zroot2
//==============================================================================
//  A lightweight value type representing an algebraic integer in the quadratic
//  extension ring  ℤ[√2]  ≔ { a + b√2 | a, b ∈ ℤ }.
//
//  Storage
//  -------
//      Integer a   coefficient of 1
//      Integer b   coefficient of √2
//
//  Example
//  -------
//      using su2compiler::ring::Zroot2;
//      Zroot2 x(3, 2);                 // 3 + 2√2
//      Zroot2 y(1, -1);                // 1 – √2
//      Zroot2 z = x * y;               // (3 + 2√2)(1 – √2) = -1 + 1√2
//      Integer norm = z.norm_sqrt();   // N(-1 + 1√2) = 1
//------------------------------------------------------------------------------
struct Zroot2
{
    Integer a = 0;
    Integer b = 0;
    
    constexpr Zroot2() noexcept = default;
    Zroot2(int _a) noexcept : a(_a), b(0) {}
    Zroot2(Integer _a) noexcept : a(_a), b(0) {}
    Zroot2(Integer _a, Integer _b) noexcept : a(_a), b(_b) {}
    
    //------------------------------------------------------------------------
    //  basic properties
    [[nodiscard]] Integer norm_sqrt2() const noexcept {
        return a * a - 2 * b * b;
    }

    [[nodiscard]] Zroot2 conj_sqrt2() const noexcept {
        return {a, -b};
    }

    template <typename T>
    [[nodiscard]] T toReal() const noexcept;

    [[nodiscard]] bool divisible(const Integer& r) const;
    [[nodiscard]] bool divisible(const Zroot2& r) const;


    //------------------------------------------------------------------------
    //  arithmetic (compound)
    Zroot2& operator+=(const Integer& r) noexcept;
    Zroot2& operator+=(const Zroot2& r) noexcept;

    Zroot2& operator-=(const Integer& r) noexcept;
    Zroot2& operator-=(const Zroot2& r) noexcept;
    Zroot2 operator-() const noexcept;

    Zroot2& operator*=(const Integer& r) noexcept;
    Zroot2& operator*=(const Zroot2& r) noexcept;

    Zroot2& operator/=(const Integer& r);
    Zroot2& operator/=(const Zroot2& r);

    auto operator<=>(const Zroot2&) const = default;
};

//------------------------------------------------------------------------
//  arithmetic (binary)
[[nodiscard]] inline Zroot2 operator+(Zroot2 lhs, const Zroot2& rhs) noexcept { return lhs += rhs; }
[[nodiscard]] inline Zroot2 operator+(Zroot2 lhs, const Integer& rhs) noexcept { return lhs += rhs; }
[[nodiscard]] inline Zroot2 operator+(const Integer& lhs, Zroot2 rhs) noexcept { return rhs += lhs; }

[[nodiscard]] inline Zroot2 operator-(Zroot2 lhs, const Zroot2& rhs) noexcept { return lhs -= rhs; }
[[nodiscard]] inline Zroot2 operator-(Zroot2 lhs, const Integer& rhs) noexcept { return lhs -= rhs; }
[[nodiscard]] inline Zroot2 operator-(const Integer& lhs, Zroot2 rhs) noexcept { return rhs -= lhs; }

[[nodiscard]] inline Zroot2 operator*(Zroot2 lhs, const Zroot2& rhs) noexcept { return lhs *= rhs; }
[[nodiscard]] inline Zroot2 operator*(Zroot2 lhs, const Integer& rhs) noexcept { return lhs *= rhs; }
[[nodiscard]] inline Zroot2 operator*(const Integer& lhs, Zroot2 rhs) noexcept { return rhs *= lhs; }

[[nodiscard]] inline Zroot2 operator/(Zroot2 lhs, const Zroot2& rhs) { return lhs /= rhs; }
[[nodiscard]] inline Zroot2 operator/(Zroot2 lhs, const Integer& rhs) { return lhs /= rhs; }

//------------------------------------------------------------------------
//  ostream
std::ostream& operator<<(std::ostream& os, const Zroot2& x);

}
}