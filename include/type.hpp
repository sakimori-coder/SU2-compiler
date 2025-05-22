#pragma once

#include <complex>
#include <iostream>
#include <Eigen/Core>
#include <gmpxx.h>
#include <mpreal.h>


namespace su2compiler
{
    using UINT = long unsigned int;

    // Define integer type
    // using Integer = long long;
    // using Integer = boost::multiprecision::int128_t;
    using Integer = mpz_class;
    using VectorXI = Eigen::Matrix<Integer, Eigen::Dynamic, 1>;
    using Vector8I = Eigen::Matrix<Integer, 8, 1>;
    using MatrixXI = Eigen::Matrix<Integer, Eigen::Dynamic, Eigen::Dynamic>;
    using Matrix8I = Eigen::Matrix<Integer, 8, 8>;

    // Define real type
    // using Real = double;
    //using Real = boost::multiprecision::cpp_dec_float_100;
    using Real = mpfr::mpreal;
    inline std::size_t PrecBits = 256;
    using VectorXR = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
    using Matrix1R = Eigen::Matrix<Real, 1, 1>;
    using Vector2R = Eigen::Matrix<Real, 2, 1>;
    using Vector4R = Eigen::Matrix<Real, 4, 1>;
    using Vector8R = Eigen::Matrix<Real, 8, 1>;
    using MatrixXR = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
    using Matrix2R = Eigen::Matrix<Real, 2, 2>; 
    using Matrix4R = Eigen::Matrix<Real, 4, 4>;
    using Matrix8R = Eigen::Matrix<Real, 8, 8>;

    // struct PrecInit {
    //     PrecInit() { mpfr::mpreal::set_default_prec(256); }
    // };
    // inline static PrecInit _mpreal_prec_init;


    // Define complex type
    using Complex = std::complex<Real>;
    using MatrixXC = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;
    using Matrix2C = Eigen::Matrix<Complex, 2, 2>;
    using Matrix4C = Eigen::Matrix<Complex, 4, 4>;


    using mpfr::sqrt;
    using mpfr::abs;
    using mpfr::log;
    using mpfr::cos;
    using mpfr::sin;

    static inline Integer round_mpreal(const Real& x) {
        Integer ret;
        mpfr_get_z(ret.get_mpz_t(), x.mpfr_ptr(), MPFR_RNDN);
        return ret;
    }

    static inline Integer ceil_mpreal(const Real& x) {
        Integer ret;
        mpfr_get_z(ret.get_mpz_t(), x.mpfr_ptr(), MPFR_RNDU);
        return ret;
    }

    static inline Integer floor_mpreal(const Real& x) {
        Integer ret;
        mpfr_get_z(ret.get_mpz_t(), x.mpfr_ptr(), MPFR_RNDD);
        return ret;
    }

    template <typename T>
    T pow_ui(T x, UINT n) {
        T ret = T(1);
        while(n > 0) {
            if(n & 1) ret *= x;
            x *= x;
            n >>= 1;
        }
        return ret;
    }



    // Constant numbers
    inline static Real PI()
    {
        return mpfr::const_pi(PrecBits);
    }

    inline static Real SQRT2()
    {
        return sqrt(Real(2, PrecBits));
    }

    inline static Real INV_SQRT2()
    {
        return Real(1, PrecBits) / SQRT2();
    }

    inline static Complex ZETA8()
    {
        return Complex(INV_SQRT2(), INV_SQRT2());
    }

    inline static Complex ZETA8_POW2()
    {
        return ZETA8() * ZETA8();
    }

    inline static Complex ZETA8_POW3()
    {
        return ZETA8() * ZETA8() * ZETA8();
    }

    inline static Complex ZETA16()
    {
        return Complex(cos(PI() / Real(8, PrecBits)), sin(PI() / Real(8, PrecBits)));
    }


    // inline std::ostream& operator<<(std::ostream& os, const mpfr::mpreal& x)
    // {
    //     std::streamsize old = os.precision();
    //     os << x.toString(os.precision());
    //     os.precision(old);
    //     return os;
    // }
}



namespace Eigen
{

    template <>
    struct NumTraits<mpfr::mpreal> : GenericNumTraits<mpfr::mpreal>
    {
        using Real        = mpfr::mpreal;
        using NonInteger  = mpfr::mpreal;
        using Nested      = mpfr::mpreal;
    
        enum {
            IsComplex             = 0,
            RequireInitialization = 1,
            ReadCost  = 8,
            AddCost   = 16,
            MulCost   = 32
        };
    

        static inline Real epsilon()
        {
            unsigned prec = mpfr::mpreal::get_default_prec();
            return mpfr::mpreal(1) >> prec;
        }
        static inline Real dummy_precision() { return Real("1e-30"); }
    
        static inline int digits10()
        {
            return int(std::floor(mpfr::mpreal::get_default_prec() * 0.3010299956639812));
        }
    };
}

namespace Eigen::numext {

    inline mpfr::mpreal sqrt (const mpfr::mpreal& x) { return mpfr::sqrt(x); }
    inline mpfr::mpreal cbrt (const mpfr::mpreal& x) { return mpfr::cbrt(x); }
    inline mpfr::mpreal exp  (const mpfr::mpreal& x) { return mpfr::exp (x); }
    inline mpfr::mpreal log  (const mpfr::mpreal& x) { return mpfr::log (x); }
    inline mpfr::mpreal sin  (const mpfr::mpreal& x) { return mpfr::sin (x); }
    inline mpfr::mpreal cos  (const mpfr::mpreal& x) { return mpfr::cos (x); }
    inline mpfr::mpreal tan  (const mpfr::mpreal& x) { return mpfr::tan (x); }
    inline mpfr::mpreal asin (const mpfr::mpreal& x) { return mpfr::asin(x); }
    inline mpfr::mpreal acos (const mpfr::mpreal& x) { return mpfr::acos(x); }
    inline mpfr::mpreal atan (const mpfr::mpreal& x) { return mpfr::atan(x); }
    inline mpfr::mpreal abs  (const mpfr::mpreal& x) { return mpfr::abs (x); }
    inline bool          isnan(const mpfr::mpreal& x){ return mpfr::isnan(x); }
    inline bool          isinf(const mpfr::mpreal& x){ return mpfr::isinf(x); }
}


namespace Eigen::internal {
    template<>
    struct cast_impl<mpz_class, mpfr::mpreal>
    {
        EIGEN_DEVICE_FUNC
        static mpfr::mpreal run(const mpz_class& z)
        {
            // mpfr::mpreal r;                              // 既定精度で確保
            // mpfr_set_z(r.mpfr_ptr(), z.get_mpz_t(),      // 最近接丸め
            //            su2compiler);
            return mpfr::mpreal(z.get_mpz_t(), su2compiler::PrecBits);                                    // NRVO でほぼゼロコピー
        }
    };
}