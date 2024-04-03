#pragma once

#include <Eigen/Core>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <complex>

namespace SU2_Compiler
{
    using ITYPE = long long;
    using UINT = long unsigned int;

    using FTYPE = boost::multiprecision::cpp_dec_float_50;
    // using FTYPE = double;

    inline static FTYPE sqrt(FTYPE x) { return boost::multiprecision::sqrt(x); }
    inline static FTYPE pow(FTYPE x, FTYPE y) { return boost::multiprecision::pow(x, y); }
    inline static FTYPE abs(FTYPE x) { return boost::multiprecision::abs(x); }

    // inline static FTYPE sqrt(FTYPE x) { return std::sqrt(x); }
    // inline static FTYPE pow(FTYPE x, FTYPE y) { return std::pow(x, y); }
    // inline static FTYPE abs(FTYPE x) { return std::abs(x); }

    using CTYPE = std::complex<FTYPE>;

    const FTYPE sqrt2 = sqrt((FTYPE)2.0);
    const FTYPE inv_sqrt2 = 1 / sqrt2;
    const CTYPE omega = {inv_sqrt2, inv_sqrt2};
    const CTYPE omega2 = omega * omega;
    const CTYPE omega3 = omega * omega * omega;
    const CTYPE sqrt_omega(cos(boost::math::constants::pi<FTYPE>() / 8), sin(boost::math::constants::pi<FTYPE>() / 8));

    using Matrix2cf = Eigen::Matrix<CTYPE, 2, 2>; 
}