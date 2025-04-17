#pragma once

#include <Eigen/Core>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <complex>

namespace SU2_Compiler
{
    using ITYPE = long long;
    // using ITYPE = __int128_t;
    // using ITYPE = __int256_t;
    // using ITYPE = boost::multiprecision::int256_t;
    using UINT = long unsigned int;

    inline static ITYPE abs(ITYPE x) { return std::abs(x); }
    // inline static ITYPE abs(ITYPE x) { return boost::multiprecision::abs(x); }

    // using FTYPE = boost::multiprecision::cpp_dec_float_50;
    using FTYPE = boost::multiprecision::number< boost::multiprecision::cpp_dec_float<50> >;
    // using FTYPE = double;

    inline static FTYPE sqrt(FTYPE x) { return boost::multiprecision::sqrt(x); }
    inline static FTYPE pow(FTYPE x, FTYPE y) { return boost::multiprecision::pow(x, y); }
    inline static FTYPE abs(FTYPE x) { return boost::multiprecision::abs(x); }
    inline static FTYPE log(FTYPE x) { return boost::multiprecision::log(x); }
    // inline static FTYPE log2(FTYPE x) { return log(x) / log((FTYPE)2.0); }
    inline static FTYPE round(FTYPE x) { return boost::multiprecision::round(x); }
    inline static FTYPE cos(FTYPE x) { return boost::multiprecision::cos(x); }
    inline static FTYPE sin(FTYPE x) { return boost::multiprecision::sin(x); }
    inline static FTYPE ceil(FTYPE x) {return boost::multiprecision::ceil(x); }
    inline static FTYPE floor(FTYPE x) {return boost::multiprecision::floor(x); }


    // inline static FTYPE sqrt(FTYPE x) { return std::sqrt(x); }
    // inline static FTYPE pow(FTYPE x, FTYPE y) { return std::pow(x, y); }
    // inline static FTYPE abs(FTYPE x) { return std::abs(x); }

    using CTYPE = std::complex<FTYPE>;

    const FTYPE PI = boost::math::constants::pi<FTYPE>();
    const FTYPE sqrt2 = sqrt((FTYPE)2.0);
    const FTYPE inv_sqrt2 = 1 / sqrt2;
    const CTYPE omega = {inv_sqrt2, inv_sqrt2};
    const CTYPE omega2 = omega * omega;
    const CTYPE omega3 = omega * omega * omega;
    const CTYPE sqrt_omega(cos(boost::math::constants::pi<FTYPE>() / 8), sin(boost::math::constants::pi<FTYPE>() / 8));
    const FTYPE LAMBDA = 1.0 + sqrt2;
    const FTYPE LAMBDA_dot = 1.0 - sqrt2;

    using Matrix2cf = Eigen::Matrix<CTYPE, 2, 2>;
    using Matrix2f = Eigen::Matrix<FTYPE, 2, 2>; 
    using Matrix4f = Eigen::Matrix<FTYPE, 4, 4>;
    using Matrix8f = Eigen::Matrix<FTYPE, 8, 8>;
    using Matrix8I = Eigen::Matrix<ITYPE, 8, 8>;
    using Vector4f = Eigen::Matrix<FTYPE, 4, 1>;
    using Vector2f = Eigen::Matrix<FTYPE, 2, 1>;
    using Vector8f = Eigen::Matrix<FTYPE, 8, 1>;
    using Vector8I = Eigen::Matrix<ITYPE, 8, 1>;

    using MatrixXf = Eigen::Matrix<FTYPE, Eigen::Dynamic, Eigen::Dynamic>;
    using MatrixXi = Eigen::Matrix<ITYPE, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXf = Eigen::Matrix<FTYPE, Eigen::Dynamic, 1>;
    using VectorXi = Eigen::Matrix<ITYPE, Eigen::Dynamic, 1>;

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