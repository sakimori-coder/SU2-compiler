#pragma once

#include <gmpxx.h>
#include <qd/dd_real.h>
#include <qd/qd_real.h>
#include <mpreal.h>
#include <complex>

#define REAL_SCALAR_TYPE_LIST \
    double,   \
    dd_real,  \
    qd_real,  \
    mpfr::mpreal

#define REAL_SCALAR_TYPES  \
    X(double)              \
    X(dd_real)             \
    X(qd_real)             \
    X(mpfr::mpreal) 


namespace su2compiler
{

using Natural = unsigned int;
using Integer = mpz_class;

template <typename T>
using Complex = std::complex<T>;


namespace type 
{

template<class... Ts> struct TypeList { };
using RealScalars = TypeList<double, dd_real, qd_real, mpfr::mpreal>;

template <typename T> T epsilon();
template <> inline double epsilon<double>() { return std::numeric_limits<double>::epsilon(); }
template <> inline dd_real epsilon<dd_real>() { return dd_real(dd_real::_eps); }
template <> inline qd_real epsilon<qd_real>() { return qd_real(qd_real::_eps); }
template <> inline mpfr::mpreal epsilon<mpfr::mpreal>() {
    unsigned p = mpfr::mpreal::get_default_prec();
    return mpfr::mpreal(1) / mpfr::pow(mpfr::mpreal(2), p);
}

template <typename T> int digits10();
template <> inline int digits10<double>() { return 15; }
template <> inline int digits10<dd_real>() { return 31; }
template <> inline int digits10<qd_real>() { return 62; }
template <> inline int digits10<mpfr::mpreal>() {
    unsigned p = mpfr::mpreal::get_default_prec();
    return static_cast<int>( std::floor(p * 0.3010299956639812) ); 
}


template <typename T> T Int_to_Real(const mpz_class& z);

template <> inline double Int_to_Real<double>(const mpz_class& z) { 
    return mpz_get_d(z.get_mpz_t());
}
template <> inline dd_real Int_to_Real<dd_real>(const mpz_class& z) {
    return dd_real(z.get_str().c_str());
}
template <> inline qd_real Int_to_Real<qd_real>(const mpz_class& z) {
    return qd_real(z.get_str().c_str());
}
template <> inline mpfr::mpreal Int_to_Real<mpfr::mpreal>(const mpz_class& z) {
    return mpfr::mpreal(z.get_mpz_t());
}


}

}

