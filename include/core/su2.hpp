#pragma once

#include <iostream>
#include <random>
#include <Eigen/Core>

#include "core/type.hpp"

namespace su2compiler
{


/**
 * @brief Structure representing a 2×2 special unitary (SU(2)) matrix.
 *
 * An SU(2) matrix is a 2×2 complex matrix U satisfying U† U = I and det(U) = 1.
 * We parametrize U by four T numbers a, b, c, and d as:
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

    SU2() noexcept : a(1), b(0), c(0), d(0) {}
    SU2(Real _a, Real _b, Real _c, Real _d) noexcept : a(_a), b(_b), c(_c), d(_d) {}
    SU2(std::complex<Real> u, std::complex<Real> t) noexcept : a(u.real()), b(u.imag()), c(t.real()), d(t.imag()) {}
    
    SU2 adjoint() const noexcept { return {a, -b, -c, -d}; }
    inline Real determinant() const noexcept { return a*a + b*b + c*c + d*d; }
    inline Real trace() const noexcept { return Real(2) * a; }
    MatrixC toMatrixC() const noexcept;
    bool isUnitary(Real tol=type::epsilon()*10) const noexcept;
    void unitalize() noexcept;

    SU2& operator*=(const SU2& r) noexcept;
};

inline SU2 operator*(SU2 lhs, const SU2& rhs) noexcept { return lhs *= rhs; }

std::ostream& operator<<(std::ostream& os, const SU2& U) noexcept;

Real distance(const SU2& U, const SU2& V) noexcept;

void set_random_unitary_seed(size_t seed) noexcept;

SU2 random_unitary() noexcept;

SU2 random_unitary(size_t seed) noexcept;

SU2 random_unitary(std::mt19937& rng) noexcept;


}