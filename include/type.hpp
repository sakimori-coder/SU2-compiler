#pragma once

#include <complex>
#include <iostream>
#include <Eigen/Core>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>


namespace su2_compiler
{
    using Integer = long long;
    // using Integer = __int128_t;
    // using Integer = __int256_t;
    // using Integer = boost::multiprecision::int256_t;
    using UINT = long unsigned int;

    // inline static Integer abs(Integer x) { return std::abs(x); }
    // inline static Integer abs(Integer x) { return boost::multiprecision::abs(x); }

    // using Real = boost::multiprecision::cpp_dec_float_50;
    using Real = boost::multiprecision::number< boost::multiprecision::cpp_dec_float<50> >;
    // using Real = double;

    inline static Real sqrt(Real x) { return boost::multiprecision::sqrt(x); }
    inline static Real pow(Real x, Real y) { return boost::multiprecision::pow(x, y); }
    inline static Real abs(Real x) { return boost::multiprecision::abs(x); }
    inline static Real log(Real x) { return boost::multiprecision::log(x); }
    // inline static Real log2(Real x) { return log(x) / log((Real)2.0); }
    inline static Real round(Real x) { return boost::multiprecision::round(x); }
    inline static Real cos(Real x) { return boost::multiprecision::cos(x); }
    inline static Real sin(Real x) { return boost::multiprecision::sin(x); }
    inline static Real ceil(Real x) {return boost::multiprecision::ceil(x); }
    inline static Real floor(Real x) {return boost::multiprecision::floor(x); }


    // inline static Real sqrt(Real x) { return std::sqrt(x); }
    // inline static Real pow(Real x, Real y) { return std::pow(x, y); }
    // inline static Real abs(Real x) { return std::abs(x); }

    using Complex = std::complex<Real>;

    const Real PI = boost::math::constants::pi<Real>();
    const Real SQRT2 = sqrt((Real)2.0);
    const Real INV_SQRT2 = 1 / SQRT2;
    const Complex ZETA8 = {INV_SQRT2, INV_SQRT2};
    const Complex ZETA8_POW2 = ZETA8 * ZETA8;
    const Complex ZETA8_POW3 = ZETA8 * ZETA8 * ZETA8;
    const Complex ZETA16(cos(boost::math::constants::pi<Real>() / 8), sin(boost::math::constants::pi<Real>() / 8));
    const Real LAMBDA = 1.0 + SQRT2;
    const Real LAMBDA_dot = 1.0 - SQRT2;

    using Matrix2C = Eigen::Matrix<Complex, 2, 2>;
    using Matrix2R = Eigen::Matrix<Real, 2, 2>; 
    using Matrix4R = Eigen::Matrix<Real, 4, 4>;
    using Matrix8R = Eigen::Matrix<Real, 8, 8>;
    using Vector4R = Eigen::Matrix<Real, 4, 1>;
    using Vector2R = Eigen::Matrix<Real, 2, 1>;
    using Vector8R = Eigen::Matrix<Real, 8, 1>;
    using Matrix8I = Eigen::Matrix<Integer, 8, 8>;
    using Vector8I = Eigen::Matrix<Integer, 8, 1>;

    using MatrixXR = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXR = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
    using MatrixXI = Eigen::Matrix<Integer, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXI = Eigen::Matrix<Integer, Eigen::Dynamic, 1>;

    inline std::ostream& operator<<( std::ostream& dest, __int128_t value )
    {
        std::ostream::sentry s( dest );
        if ( s ) {
            __uint128_t tmp = value < 0 ? -value : value;
            char buffer[ 128 ];
            char* d = std::end( buffer );
            do
            {
                -- d;
                *d = "0123456789"[ tmp % 10 ];
                tmp /= 10;
            } while ( tmp != 0 );
            if ( value < 0 ) {
                -- d;
                *d = '-';
            }
            int len = std::end( buffer ) - d;
            if ( dest.rdbuf()->sputn( d, len ) != len ) {
                dest.setstate( std::ios_base::badbit );
            }
        }
        return dest;
    }
}