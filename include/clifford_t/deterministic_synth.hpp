#pragma once

#include <string>
#include <vector>

#include "core/type.hpp"
#include "core/su2.hpp"
#include "ring/Zzeta8j.hpp"
#include "clifford_t/exact_synth.hpp"

namespace su2compiler::clifford_t::deterministic {

std::string synth(
        SU2 V,
        Real eps
);

std::vector<exact::U2Dzeta8> fixed_t_synth(
        SU2 V,
        Real eps,
        Natural t
);


std::vector<ring::Zzeta8j> solve_approx_lattice(
        SU2 V,
        Real eps,
        Natural k
);

}