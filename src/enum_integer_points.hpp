#pragma once

#include "type.hpp"

namespace SU2_Compiler{

/**
 * @brief Performs the LLL (Lenstra–Lenstra–Lovász) lattice basis reduction.
 *
 * Given an input basis matrix B whose columns are the original lattice basis vectors,
 * this function returns a reduced integer basis matrix that satisfies the LLL
 * reduction conditions under the parameter δ.
 *
 * @param[in] B      Floating-point m×n basis matrix
 * @param[in] delta  Reduction parameter (0.25 < δ < 1.0)
 * @return           n×n unimodular integer matrix U such that B * U is an LLL-reduced basis.
 *
 * @see
 *     https://en.wikipedia.org/wiki/LLL_algorithm
 */
MatrixXi LLL(MatrixXf B, FTYPE delta);


ITYPE determinant(MatrixXi A);
MatrixXi inverse_Unimodular(MatrixXi U);
std::vector<VectorXi> EnumIntegerPoints(MatrixXf Q, VectorXf p, FTYPE c);
    

}