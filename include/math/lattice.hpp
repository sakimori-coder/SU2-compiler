#pragma once

#include <utility>
#include <Eigen/Core>

#include "core/type.hpp"

namespace su2compiler::math::lattice{

/**
 * @brief Performs the LLL (Lenstra–Lenstra–Lovász) lattice basis reduction.
 *
 * Given an input basis matrix B whose columns are the original lattice basis vectors,
 * this function returns a reduced integer basis matrix that satisfies the LLL
 * reduction conditions under the parameter δ.
 *
 * @param[in] B      Floating-point m×n basis matrix
 * @param[in] delta  Reduction parameter (0.25 < δ < 1.0)
 * @return           A pair of n×n integer matrices (U, U_inv), where U is a Unimodular matrix 
 *                   such that B * U is an LLL-reduced basis and its inverse.
 *
 * @see
 *     https://en.wikipedia.org/wiki/LLL_algorithm
 */
std::pair<MatrixI, MatrixI> LLL(
        MatrixR B, 
        Real delta
);


/**
 * @brief Enumerate all integer points within the ellipsoid defined by (x - p)^T Q (x - p) <= c.
 *
 * This function finds and returns all integer lattice points x ∈ ℤ^n that lie inside or on the boundary
 * of the ellipsoid centered at vector p with shape matrix Q and radius squared c.
 *
 * The ellipsoid is defined by the inequality:
 *     (x - p)^T Q (x - p) <= c
 * where Q is a symmetric positive definite matrix.
 *
 * @param Q    A positive definite real symmetric matrix defining the ellipsoid shape.
 * @param p    A real vector representing the center of the ellipsoid.
 * @param c    A real scalar representing the squared "radius" of the ellipsoid.
 * @return     A vector of integer vectors (VectorXI) containing all points satisfying the inequality.
 */
std::vector<VectorI> EnumIntegerPoints(
        MatrixR Q,
        VectorR p,
        Real c
);
    

}