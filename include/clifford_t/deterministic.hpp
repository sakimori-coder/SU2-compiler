#pragma once

#include <string>
#include <vector>

#include "ring/Zzeta8j.hpp"
#include "type.hpp"
#include "su2.hpp"
#include "clifford_t/exact_synthesis.hpp"

namespace su2compiler {
namespace clifford_t {


std::string deterministic_synthesis(SU2 V, Real eps);
std::vector<U2Dzeta8> deterministic_synthesis_fixed_t(SU2 V, Real eps, UINT t);
std::vector<ring::Zzeta8j> solve_approx_synthesis(SU2 V, Real eps, UINT k);


}
}