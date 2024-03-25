#pragma once

#include <string>
#include <vector>
#include <Eigen/Core>
#include "type.hpp"
#include "rings.hpp"
#include "quaternion.hpp"
#include "U2_ZOmega.hpp"

namespace SU2_Compiler
{
    using Mat3ZRoot2 = Eigen::Matrix<ZRoot2, 3, 3>;

    std::pair<Mat3ZRoot2, int> U2_to_SO3(const U2_ZOmega& U);
    Eigen::Matrix3i parity(const Mat3ZRoot2& A);
    std::string ExactSynthesis(const U2_ZOmega& U);
    int get_T_count(const U2_ZOmega& U);
    quaternion to_quaternion(std::string);
}