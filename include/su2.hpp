#pragma once

#include <iostream>
#include <random>
#include <Eigen/Core>

#include "type.hpp"

namespace su2compiler
{


/**
 * @brief Structure representing a 2×2 special unitary (SU(2)) matrix.
 *
 * An SU(2) matrix is a 2×2 complex matrix U satisfying U† U = I and det(U) = 1.
 * We parametrize U by four real numbers a, b, c, and d as:
 *
 *   U = [ a + i b   ,  –c + i d ]
 *       [ c + i d   ,   a – i b ]
 *
 * This parametrization automatically enforces the normalization
 *   a² + b² + c² + d² = 1.
 */
struct SU2
{
    Real a;
    Real b;
    Real c;
    Real d;

    SU2() : a(1), b(0), c(0), d(0) {}
    SU2(Real _a, Real _b, Real _c, Real _d) : a(_a), b(_b), c(_c), d(_d) {}
    SU2(Complex u, Complex t) : a(u.real()), b(u.imag()), c(t.real()), d(t.imag()) {}
    
    SU2 adjoint() const { return {a, -b, -c, -d}; }
    inline Real determinant() const { return a*a + b*b + c*c + d*d; }
    inline Real trace() const { return 2*a; }
    Matrix2C to_EigenMatrix() const;
    bool is_unitary(Real to=1e-14) const;
    void unitalize();

    SU2& operator*=(const SU2& r);
    SU2 operator*(const SU2& r) const;
};

std::ostream& operator<<(std::ostream& os, const SU2& U);

Real distance(const SU2& U, const SU2& V);
void set_random_unitary_seed(size_t seed);
SU2 random_unitary();
SU2 random_unitary(size_t seed);
SU2 random_unitary(std::mt19937& rng);


}