#pragma once

#include "type.hpp"
#include "su2.hpp"


namespace su2compiler::mixing_su2
{


MatrixXC Choi_Jamiolkowski(const MatrixXC& U);
Matrix4R Choi_Jamiolkowski_MagicBasis(const Matrix2C& U);
std::pair<Real, std::vector<Real>> optimize_distribution(const std::vector<SU2>& availableU, SU2 V);


}