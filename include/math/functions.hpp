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
inline double       sqrt(double  x) noexcept { return std::sqrt(x); }
inline dd_real      sqrt(dd_real x) noexcept { return    ::sqrt(x); }
inline qd_real      sqrt(qd_real x) noexcept { return    ::sqrt(x); }
inline mpfr::mpreal sqrt(const mpfr::mpreal&  x) noexcept { return mpfr::sqrt(x); } 

// abs
inline double       abs(double  x) noexcept { return std::abs(x); }
inline dd_real      abs(dd_real x) noexcept { return    ::abs(x); }
inline qd_real      abs(qd_real x) noexcept { return    ::abs(x); }
inline mpfr::mpreal abs(const mpfr::mpreal& x) noexcept { return mpfr::abs(x); }
template <typename T>
inline T abs(std::complex<T> x) noexcept {
     return math::sqrt(x.real() * x.real() + x.imag() * x.imag()); 
}

// log
inline double       log(double  x) noexcept { return std::log(x); }
inline dd_real      log(dd_real x) noexcept { return    ::log(x); }
inline qd_real      log(qd_real x) noexcept { return    ::log(x); }
inline mpfr::mpreal log(const mpfr::mpreal& x) noexcept { return mpfr::log(x); }

// log2
inline double       log2(double  x) noexcept { return math::log(x) / math::log(2.0); }
inline dd_real      log2(dd_real x) noexcept { return math::log(x) / math::log(dd_real(2)); }
inline qd_real      log2(qd_real x) noexcept { return math::log(x) / math::log(qd_real(2)); }
inline mpfr::mpreal log2(const mpfr::mpreal& x) noexcept { return math::log(x) / math::log(mpfr::mpreal(2)); }

// cos
inline double       cos(double  x) noexcept { return std::cos(x); }
inline dd_real      cos(dd_real x) noexcept { return    ::cos(x); }
inline qd_real      cos(qd_real x) noexcept { return    ::cos(x); }
inline mpfr::mpreal cos(const mpfr::mpreal& x) noexcept { return mpfr::cos(x); }

// sin
inline double       sin(double  x) noexcept { return std::sin(x); }
inline dd_real      sin(dd_real x) noexcept { return    ::sin(x); }
inline qd_real      sin(qd_real x) noexcept { return    ::sin(x); }
inline mpfr::mpreal sin(const mpfr::mpreal& x) noexcept { return mpfr::sin(x); }

// round (output integer)
inline mpz_class round(double  x) noexcept { return mpz_class(std::round(x)); }
inline mpz_class round(dd_real x) noexcept { return math::round(to_double(x)); }
inline mpz_class round(qd_real x) noexcept { return math::round(to_double(x)); }
inline mpz_class round(const mpfr::mpreal& x) noexcept {
    mpz_class ret;
    mpfr_get_z(ret.get_mpz_t(), x.mpfr_ptr(), MPFR_RNDN);
    return ret;
}

// ceil
inline mpz_class ceil(double  x) noexcept { return mpz_class(std::ceil(x)); }
inline mpz_class ceil(dd_real x) noexcept { return math::ceil(to_double(x)); }
inline mpz_class ceil(qd_real x) noexcept { return math::ceil(to_double(x)); }
inline mpz_class ceil(const mpfr::mpreal& x) noexcept {
    mpz_class ret;
    mpfr_get_z(ret.get_mpz_t(), x.mpfr_ptr(), MPFR_RNDU);
    return ret;
}

// floor
inline mpz_class floor(double  x) noexcept { return mpz_class(std::floor(x)); }
inline mpz_class floor(dd_real x) noexcept { return math::floor(to_double(x)); }
inline mpz_class floor(qd_real x) noexcept { return math::floor(to_double(x)); }
inline mpz_class floor(const mpfr::mpreal& x) noexcept {
    mpz_class ret;
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
[[gnu::always_inline]] inline void FMA(double x, double y, double z, double& res) noexcept { 
    res = std::fma(x,y,z);
}
[[gnu::always_inline]] inline void FMA(dd_real x, dd_real y, dd_real z, dd_real& res) noexcept {
    res = x * y + z;
}
[[gnu::always_inline]] inline void FMA(qd_real x, qd_real y, qd_real z, qd_real& res) noexcept {
    res = x * y + z;
}
[[gnu::always_inline]] inline void FMA(
    const mpfr::mpreal& x, 
    const mpfr::mpreal& y, 
    const mpfr::mpreal& z,
    mpfr::mpreal& res
) noexcept 
{
    mpfr_fma(res.mpfr_ptr(), x.mpfr_ptr(), y.mpfr_ptr(), z.mpfr_ptr(), MPFR_RNDN);
}

}