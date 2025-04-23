#include "su2.hpp"

#include <iostream>
#include <complex>
#include <random>
#include <Eigen/Core>

#include "type.hpp"

namespace su2_compiler
{


Matrix2C SU2::to_EigenMatrix() const
{
    Matrix2C ret;
    Complex u(a,b);
    Complex t(c,d);
    ret << u, -std::conj(t),
           t,  std::conj(u);
    return ret;
}


bool SU2::is_unitary(Real tol) const
{
    if(abs(this->determinant() - (Real)1.0) < tol) return true;
    else return false; 
}


void SU2::unitalize()
{
    Real norm = sqrt(this->determinant());
    a /= norm;
    b /= norm;
    c /= norm;
    d /= norm;
}


SU2& SU2::operator*=(const SU2& r)
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


SU2 SU2::operator*(const SU2& r) const
{
    SU2 ret = (*this);
    ret *= r;
    return ret;
}


std::ostream& operator<<(std::ostream& os, const SU2& U)
{
    os << "(" << U.a << ", " << U.b << ", " << U.c << ", " << U.d << ")";
    return os;
}


// SU2 adjoint(const SU2& U)
// {
//     return {U.a, -U.b, -U.c, -U.d};
// }


Real distance(const SU2& U, const SU2& V)
{
    SU2 UV_dag = U * V.adjoint();
    Real tr = UV_dag.trace();
    return sqrt(1.0 - (tr*tr) / 4.0);
}


SU2 random_unitary()
{
    static std::mt19937 rng{ std::random_device{}() };
    return random_unitary(rng);
}


SU2 random_unitary(size_t seed)
{
    std::mt19937 rng(seed);
    return random_unitary(rng);
}


SU2 random_unitary(std::mt19937& rng)
{
    std::normal_distribution<double> gauss(0.0, 1.0);
    Real a = (Real)gauss(rng);
    Real b = (Real)gauss(rng);
    Real c = (Real)gauss(rng);
    Real d = (Real)gauss(rng);

    SU2 ret(a,b,c,d);
    ret.unitalize();
    return ret;
}

}