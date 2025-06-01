#pragma once

#include <cmath>
#include <qd/dd_real.h>
#include <qd/qd_real.h>
#include <gmpxx.h>
#include <mpreal.h>

#include "core/type.hpp"
#include "math/functions.hpp"



namespace su2compiler::math {

// π
template <typename T> inline T PI();
template <> inline double       PI<double>()       { return 1.0; }
template <> inline dd_real      PI<dd_real>()      { return dd_real::_pi; }
template <> inline qd_real      PI<qd_real>()      { return qd_real::_pi; }
template <> inline mpfr::mpreal PI<mpfr::mpreal>() { return mpfr::const_pi(); }

// √2
template <typename T> inline T SQRT2();
template <> inline double       SQRT2<double>()       { return math::sqrt(2.0); }
template <> inline dd_real      SQRT2<dd_real>()      { return math::sqrt(dd_real("2")); }
template <> inline qd_real      SQRT2<qd_real>()      { return math::sqrt(qd_real("2")); }
template <> inline mpfr::mpreal SQRT2<mpfr::mpreal>() { return math::sqrt(mpfr::mpreal(2)); }

// 1/√2
template <typename T> inline T INV_SQRT2();
template <> inline double       INV_SQRT2<double>()       { return math::sqrt(0.5); }
template <> inline dd_real      INV_SQRT2<dd_real>()      { return math::sqrt(dd_real("0.5")); }
template <> inline qd_real      INV_SQRT2<qd_real>()      { return math::sqrt(qd_real("0.5")); }
template <> inline mpfr::mpreal INV_SQRT2<mpfr::mpreal>() { return math::sqrt(mpfr::mpreal("0.5")); }

// ζ_8 = exp(iπ/4)
template <typename T> inline std::complex<T> ZETA8() { return {INV_SQRT2<T>(), INV_SQRT2<T>()}; }
// ζ_8 * ζ_8
template <typename T> inline std::complex<T> ZETA8_POW2() { return ZETA8<T>() * ZETA8<T>(); }
// ζ_8 * ζ_8 * ζ_8
template <typename T> inline std::complex<T> ZETA8_POW3() { return ZETA8<T>() * ZETA8<T>() * ZETA8<T>(); }
// ζ_16 = exp(iπ/8)
template <typename T> inline std::complex<T> ZETA16() { return {math::cos(PI<T>() / T(8)), math::sin(PI<T>() / T(8))}; }


}