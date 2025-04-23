#pragma once

#include <string>
#include <Eigen/Core>

#include "include/type.hpp"
#include "include/rings.hpp"
#include "include/SpecialUnitary2.hpp"

namespace SU2_Compiler
{
    using Matrix3ZRoot2 = Eigen::Matrix<ZRoot2, 3, 3>;

    std::string ExactSynthesis(const DOmega_U2& U);
    int get_T_count(const DOmega_U2& U);
    SpecialUnitary2 String_to_SU2(std::string);
}