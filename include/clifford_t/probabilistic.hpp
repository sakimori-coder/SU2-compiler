#pragma once

#include <string>
#include <utility>
#include <vector>

#include "type.hpp"
#include "su2.hpp"

namespace su2compiler::clifford_t
{


std::vector<std::pair<Real, std::string>> probabilistic_synthesis(SU2 V, Real eps);
std::vector<std::pair<Real, std::string>> probabilistic_synthesis_fixed_t(SU2 V, Real eps, UINT t);


}