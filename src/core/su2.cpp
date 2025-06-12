#include "core/su2.hpp"

#include <iostream>
#include <random>
#include <Eigen/Core>

#include "core/type.hpp"
#include "math/functions.hpp"


namespace su2compiler
{


MatrixC SU2::toMatrixC() const noexcept
{
    MatrixC ret(2,2);
    Complex u(a,b);
    Complex t(c,d);
    ret << u, -std::conj(t),
           t,  std::conj(u);
    return ret;
}


bool SU2::isUnitary(Real tol) const noexcept
{
    if(math::abs(this->determinant() - Real(1)) < type::epsilon()*10) return true;
    else return false; 
}


void SU2::unitalize() noexcept
{
    Real norm = math::sqrt(this->determinant());
    a /= norm;
    b /= norm;
    c /= norm;
    d /= norm;
}


SU2& SU2::operator*=(const SU2& r) noexcept
{
        Real ret_a = a*r.a - b*r.b - c*r.c - d*r.d;
        Real ret_b = a*r.b + b*r.a - c*r.d + d*r.c;
        Real ret_c = a*r.c + b*r.d + c*r.a - d*r.b;
        Real ret_d = a*r.d - b*r.c + c*r.b + d*r.a;
    
        a = ret_a;
        b = ret_b;
        c = ret_c;
        d = ret_d;

        return *this;
}


std::ostream& operator<<(std::ostream& os, const SU2& U) noexcept
{
    os << "(" << U.a << ", " << U.b << ", " << U.c << ", " << U.d << ")";
    return os;
}


Real distance(const SU2& U, const SU2& V) noexcept
{
    SU2 UV_dag = U * V.adjoint();
    Real tr = UV_dag.trace();
    return math::sqrt(Real(1) - (tr*tr) / Real(4));
}


SU2 random_unitary() noexcept
{
    static std::mt19937 rng{ std::random_device{}() };
    return random_unitary(rng);
}


SU2 random_unitary(size_t seed) noexcept
{
    std::mt19937 rng(seed);
    return random_unitary(rng);
}


SU2 random_unitary(std::mt19937& rng) noexcept
{
    std::normal_distribution<double> gauss(0.0, 1.0);
    Real a = Real(gauss(rng));
    Real b = Real(gauss(rng));
    Real c = Real(gauss(rng));
    Real d = Real(gauss(rng));

    SU2 ret(a,b,c,d);
    ret.unitalize();
    return ret;
}


}