#pragma once

#include <string>
#include <vector>

#include "ring/Zzeta8j.hpp"
#include "type.hpp"
#include "su2.hpp"
#include "clifford_t/exact_synthesize.hpp"

namespace su2_compiler {
namespace clifford_t {


std::string synthesize(SU2 V, Real eps);
std::vector<ring::Zzeta8j> solve_approx_synthesis(SU2 V, Real eps, int k);


}
}