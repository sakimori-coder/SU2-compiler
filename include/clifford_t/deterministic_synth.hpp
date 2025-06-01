#pragma once

#include <string>
#include <vector>

#include "core/type.hpp"
#include "core/su2.hpp"
#include "ring/Zzeta8j.hpp"
#include "clifford_t/exact_synth.hpp"

namespace su2compiler::clifford_t::deterministic {

template <typename RealType>
std::string synth(
        SU2<RealType> V,
        RealType eps
);

template <typename RealType>
std::vector<exact::U2Dzeta8> fixed_t_synth(
        SU2<RealType> V,
        RealType eps,
        Natural t
);

template <typename RealType>
std::vector<ring::Zzeta8j> solve_approx_lattice(
        SU2<RealType> V,
        RealType eps,
        Natural k
);

}