#include "core/su2.hpp"

#include <iostream>
#include <random>
#include <Eigen/Core>

#include "core/type.hpp"
#include "math/functions.hpp"


namespace su2compiler
{

template <typename RealType>
Eigen::Matrix2<std::complex<RealType>> SU2<RealType>::toEigenMatrix() const noexcept
{
    Eigen::Matrix2<std::complex<RealType>> ret;
    std::complex<RealType> u(a,b);
    std::complex<RealType> t(c,d);
    ret << u, -std::conj(t),
           t,  std::conj(u);
    return ret;
}
#define X(RealType) template Eigen::Matrix2<std::complex<RealType>> SU2<RealType>::toEigenMatrix() const noexcept;
REAL_SCALAR_TYPES
#undef X


template <typename T>
bool SU2<T>::isUnitary(T tol) const noexcept
{
    if(math::abs(this->determinant() - T(1)) < type::epsilon<T>()*10) return true;
    else return false; 
}
#define X(T) template bool SU2<T>::isUnitary(T) const noexcept;
REAL_SCALAR_TYPES
#undef X


template <typename T>
void SU2<T>::unitalize() noexcept
{
    T norm = math::sqrt(this->determinant());
    a /= norm;
    b /= norm;
    c /= norm;
    d /= norm;
}
#define X(T) template void SU2<T>::unitalize() noexcept;
REAL_SCALAR_TYPES
#undef X


template <typename T>
SU2<T>& SU2<T>::operator*=(const SU2<T>& r) noexcept
{
        T ret_a = a*r.a - b*r.b - c*r.c - d*r.d;
        T ret_b = a*r.b + b*r.a - c*r.d + d*r.c;
        T ret_c = a*r.c + b*r.d + c*r.a - d*r.b;
        T ret_d = a*r.d - b*r.c + c*r.b + d*r.a;
    
        a = ret_a;
        b = ret_b;
        c = ret_c;
        d = ret_d;

        return *this;
}
#define X(T) template SU2<T>& SU2<T>::operator*=(const SU2<T>&) noexcept;
REAL_SCALAR_TYPES
#undef X


template <typename T>
std::ostream& operator<<(std::ostream& os, const SU2<T>& U) noexcept
{
    os << "(" << U.a << ", " << U.b << ", " << U.c << ", " << U.d << ")";
    return os;
}
#define X(T) template std::ostream& operator<<(std::ostream&, const SU2<T>&) noexcept;
REAL_SCALAR_TYPES
#undef X


template <typename T>
T distance(const SU2<T>& U, const SU2<T>& V) noexcept
{
    SU2<T> UV_dag = U * V.adjoint();
    T tr = UV_dag.trace();
    return math::sqrt(T(1) - (tr*tr) / T(4));
}
#define X(T) template T distance(const SU2<T>&, const SU2<T>&) noexcept;
REAL_SCALAR_TYPES
#undef X


template <typename T>
SU2<T> random_unitary() noexcept
{
    static std::mt19937 rng{ std::random_device{}() };
    return random_unitary<T>(rng);
}
#define X(T) template SU2<T> random_unitary() noexcept;
REAL_SCALAR_TYPES
#undef X


template <typename T>
SU2<T> random_unitary(size_t seed) noexcept
{
    std::mt19937 rng(seed);
    return random_unitary<T>(rng);
}
#define X(T) template SU2<T> random_unitary(size_t) noexcept;
REAL_SCALAR_TYPES
#undef X


template <typename T>
SU2<T> random_unitary(std::mt19937& rng) noexcept
{
    std::normal_distribution<double> gauss(0.0, 1.0);
    T a = T(gauss(rng));
    T b = T(gauss(rng));
    T c = T(gauss(rng));
    T d = T(gauss(rng));

    SU2<T> ret(a,b,c,d);
    ret.unitalize();
    return ret;
}
#define X(T) template SU2<T> random_unitary(std::mt19937&) noexcept;
REAL_SCALAR_TYPES
#undef X


}