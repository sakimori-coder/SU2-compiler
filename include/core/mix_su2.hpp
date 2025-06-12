#pragma once

#include <vector>
#include <Eigen/Core>

#include "core/type.hpp"
#include "core/su2.hpp"
#include "math/sdp.hpp"

namespace su2compiler
{

class MixSU2
{
public:
    VectorR prob;
    std::vector<SU2> availableU;

    MixSU2() = default;
    MixSU2(const Eigen::VectorX<Real>& _prob, 
           const std::vector<SU2>& _availableU) noexcept : prob(_prob), availableU(_availableU) {}
    MixSU2(const std::vector<SU2>& _availableU) noexcept : availableU(_availableU) {}

    math::sdp::Results compute_optimal_prob(
            const SU2& targetV,
            math::sdp::Options opt = math::sdp::Options{}
    );

};



}