#pragma once

#include <string>
#include <Eigen/Core>

#include "type.hpp"
#include "ring/Zroot2.hpp"
#include "ring/Zzeta8j.hpp"

namespace su2_compiler {
namespace clifford_t {

using Matrix2Zzeta8 = Eigen::Matrix<ring::Zzeta8,2,2>;
using Matrix3Zroot2 = Eigen::Matrix<ring::Zroot2,3,3>;

std::string exact_synthesis(ring::Zzeta8 u, ring::Zzeta8 t, int l, Integer k);
std::string exact_synthesis(const std::string& sequence);

}
}