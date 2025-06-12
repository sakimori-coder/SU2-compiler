#pragma once

#include <cmath>

#include "core/type.hpp"
#include "math/functions.hpp"



namespace su2compiler::math {

// π
inline Real PI() { return mpfr::const_pi(); }

// √2
inline Real SQRT2() { return math::sqrt(Real(2)); }

// 1/√2
inline Real INV_SQRT2() { return math::sqrt(Real("0.5")); }

// ζ_8 = exp(iπ/4)
inline Complex ZETA8() { return {INV_SQRT2(), INV_SQRT2()}; }
// ζ_8 * ζ_8
inline Complex ZETA8_POW2() { return ZETA8() * ZETA8(); }
// ζ_8 * ζ_8 * ζ_8
inline Complex ZETA8_POW3() { return ZETA8() * ZETA8() * ZETA8(); }
// ζ_16 = exp(iπ/8)
inline Complex ZETA16() { return {math::cos(PI() / Real(8)), math::sin(PI() / Real(8))}; }


}