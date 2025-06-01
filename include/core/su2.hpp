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
template <typename RealType>
struct SU2
{
    RealType a;
    RealType b;
    RealType c;
    RealType d;

    SU2() noexcept : a(1), b(0), c(0), d(0) {}
    SU2(RealType _a, RealType _b, RealType _c, RealType _d) noexcept : a(_a), b(_b), c(_c), d(_d) {}
    SU2(std::complex<RealType> u, std::complex<RealType> t) noexcept : a(u.real()), b(u.imag()), c(t.real()), d(t.imag()) {}
    
    SU2<RealType> adjoint() const noexcept { return {a, -b, -c, -d}; }
    inline RealType determinant() const noexcept { return a*a + b*b + c*c + d*d; }
    inline RealType trace() const noexcept { return RealType(2) * a; }
    Eigen::Matrix2<std::complex<RealType>> toEigenMatrix() const noexcept;
    bool isUnitary(RealType tol=type::epsilon<RealType>()*10) const noexcept;
    void unitalize() noexcept;

    SU2& operator*=(const SU2& r) noexcept;
};


template <typename T>
inline SU2<T> operator*(SU2<T> lhs, const SU2<T>& rhs) noexcept { return lhs *= rhs; }


template <typename T>
std::ostream& operator<<(std::ostream& os, const SU2<T>& U) noexcept;

template <typename T>
T distance(const SU2<T>& U, const SU2<T>& V) noexcept;

void set_random_unitary_seed(size_t seed) noexcept;

template <typename T>
SU2<T> random_unitary() noexcept;

template <typename T>
SU2<T> random_unitary(size_t seed) noexcept;

template <typename T>
SU2<T> random_unitary(std::mt19937& rng) noexcept;


}