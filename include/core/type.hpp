#pragma once

#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <mpreal.h>
#include <complex>
#include <Eigen/Core>


namespace su2compiler
{

using Real = mpfr::mpreal;
using Complex = std::complex<Real>;
using Natural = unsigned int;
using Integer = mpz_class;

using MatrixR = Eigen::MatrixX<Real>;
using VectorR = Eigen::VectorX<Real>;

using MatrixC = Eigen::MatrixX<Complex>;
using VectorC = Eigen::VectorX<Complex>;

using MatrixI = Eigen::MatrixX<Integer>;
using VectorI = Eigen::VectorX<Integer>;


enum class STATUS{
    SUCCESS,
    FAILURE,
    TIMEOUT,
};


namespace type 
{


inline mpfr::mpreal epsilon() {
    unsigned p = mpfr::mpreal::get_default_prec();
    return mpfr::mpreal(1) / mpfr::pow(mpfr::mpreal(2), p);
}

inline int digits10() {
    unsigned p = mpfr::mpreal::get_default_prec();
    return static_cast<int>( std::floor(p * 0.3010299956639812) ); 
}



inline mpfr::mpreal Int_to_Real(const mpz_class& z) {
    return mpfr::mpreal(z.get_mpz_t());
}


}

}

