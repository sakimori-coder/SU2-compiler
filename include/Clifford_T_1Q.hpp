#pragma once

#include <string>
#include <Eigen/Core>

#include "type.hpp"
#include "rings.hpp"
#include "SpecialUnitary2.hpp"

namespace SU2_Compiler
{

using Matrix2ZOmega = Eigen::Matrix<ZOmega, 2, 2>;
using Matrix3ZRoot2 = Eigen::Matrix<ZRoot2, 3, 3>;

extern ZOmega omega_pow[8];

class Clifford_T_1Q
{
private:
    Matrix2ZOmega mat;
    Integer k;

public:
    Clifford_T_1Q();
    Clifford_T_1Q(Eigen::Matrix<ZOmega, 2, 2> _mat, Integer _k) : mat(_mat), k(_k) {};
    Clifford_T_1Q(const std::string& sequence);
    Clifford_T_1Q(ZOmega u, ZOmega t, int l, Integer k);

    std::string compute_sequence();
    int compute_Tcount() const;
    SpecialUnitary2 to_SU2() const;
};


}
