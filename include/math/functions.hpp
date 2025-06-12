#pragma once

#include <cmath>
#include <cstdint>
#include <qd/dd_real.h>
#include <qd/qd_real.h>
#include <gmpxx.h>
#include <mpreal.h>

#include "core/type.hpp"

namespace su2compiler::math {

// sqrt
inline Real sqrt(const Real&  x) noexcept { return mpfr::sqrt(x); } 

// abs
inline Real abs(const Real& x) noexcept { return mpfr::abs(x); }
inline Real abs(const Complex& x) noexcept {
    return math::sqrt(x.real() * x.real() + x.imag() * x.imag());
}

// log
inline Real log(const Real& x) noexcept { return mpfr::log(x); }

// log2
inline Real log2(const Real& x) noexcept { return math::log(x) / math::log(Real(2)); }

// cos
inline Real cos(const Real& x) noexcept { return mpfr::cos(x); }

// sin
inline Real sin(const Real& x) noexcept { return mpfr::sin(x); }

// round (output integer)
inline Integer round(const Real& x) noexcept {
    Integer ret;
    mpfr_get_z(ret.get_mpz_t(), x.mpfr_ptr(), MPFR_RNDN);
    return ret;
}

// ceil
inline Integer ceil(const Real& x) noexcept {
    Integer ret;
    mpfr_get_z(ret.get_mpz_t(), x.mpfr_ptr(), MPFR_RNDU);
    return ret;
}

// floor
inline Integer floor(const Real& x) noexcept {
    Integer ret;
    mpfr_get_z(ret.get_mpz_t(), x.mpfr_ptr(), MPFR_RNDD);
    return ret;
}

// pow (x^n)
template <typename T>
T pow_ui(T x, std::uint32_t n) noexcept {
    T ret = T(1);
    while(n > 0) {
        if(n & 1) ret *= x;
        x *= x;
        n >>= 1;
    }
    return ret;
}

// FMA (Fused Multiply-Add) 
// res = x * y + z
[[gnu::always_inline]] inline void FMA(
    Real& res,
    const Real& x, 
    const Real& y,
    const Real& z
) noexcept 
{
    mpfr_fma(res.mpfr_ptr(), x.mpfr_ptr(), y.mpfr_ptr(), z.mpfr_ptr(), MPFR_RNDN);
}

}