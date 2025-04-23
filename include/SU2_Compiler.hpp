#pragma once

#include <vector>
#include <tuple>
#include <string>

#include "type.hpp"
#include "rings.hpp"
#include "SpecialUnitary2.hpp"

namespace SU2_Compiler
{

std::string Deterministic_Clifford_T(SpecialUnitary2 V, Real eps);

std::vector<std::tuple<ZOmega, ZOmega, Integer>> solve_approximation_synthesis_CliffordT(SpecialUnitary2 V, Real eps);


}