#pragma once

#include <tuple>

#include "type.hpp"


namespace su2_compiler{



/**
 * @brief Performs Gram–Schmidt orthogonalization on the columns of a basis matrix.
 *
 * Given an input basis matrix \p B whose columns are the original (possibly non-orthogonal)
 * basis vectors, this function returns a matrix of the same dimensions whose columns
 * form an orthonormal basis computed via the Gram–Schmidt process.
 *
 * @param[in] B  
 *     An m×n matrix whose n columns are the input basis vectors in ℝᵐ.
 * @return 
 *     An m×n matrix whose columns are the orthogonalized vectors corresponding to \p B.
 *     Note: columns are not scaled to unit length.
 *
 * @see  
 *     https://en.wikipedia.org/wiki/Gram–Schmidt_process
 */
MatrixXR GSO(const MatrixXR& B);



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
std::tuple<MatrixXI, MatrixXI> LLL(MatrixXR B, Real delta);



/**
 * @brief Compute the Cholesky decomposition of a symmetric positive definite matrix A.
 *
 * This function returns the lower-triangular matrix L such that:
 *     A = L * L.transpose()
 * where A is a real symmetric positive definite matrix.
 *
 * @param A    A symmetric positive definite matrix.
 * @return     A lower-triangular matrix L such that A = L * L^T.
 */
MatrixXR cholesky(const MatrixXR& A);



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
std::vector<VectorXI> EnumIntegerPoints(MatrixXR Q, VectorXR p, Real c);
    

}