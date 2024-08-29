#pragma once

#include <vector>
#include <utility>
#include <unordered_map>
#include <tbb/concurrent_hash_map.h>
#include "type.hpp"
#include "quaternion.hpp"
#include "U2_ZOmega.hpp"
#include "HashTable.hpp"


namespace SU2_Compiler
{
    std::pair< std::vector< U2_ZOmega >, std::array<std::vector< ZRoot2 >, 12> >
    enum_u_t(quaternion U, FTYPE eps, int k, int l, const std::array<std::vector< ZRoot2 >, 12>& pre_results = {});

    std::vector<U2_ZOmega>
    enum_u_t
    (quaternion U, FTYPE eps, FTYPE eps_pre, int k, int l, std::vector<ZRoot2>& XY, std::vector<ZRoot2>& ZW, std::vector<U2_ZOmega>& candidates);

    std::vector<U2_ZOmega> enum_u_t_use_hash
    (quaternion U, FTYPE eps, FTYPE eps_pre, int k, int l, tbb::concurrent_hash_map<ZRoot2, std::pair<ZRoot2, ZRoot2>, ZRoot2_hash>& hash_XY, tbb::concurrent_hash_map<ZRoot2, std::pair<ZRoot2, ZRoot2>, ZRoot2_hash>& hash_ZW);
}