#pragma once

#include <vector>
#include <utility>
#include "type.hpp"
#include "quaternion.hpp"
#include "U2_ZOmega.hpp"


namespace SU2_Compiler
{
    std::pair< std::vector< U2_ZOmega >, std::array<std::vector< ZRoot2 >, 12> >
    enum_u_t(quaternion U, FTYPE eps, int k, int l, const std::array<std::vector< ZRoot2 >, 12>& pre_results = {});
}