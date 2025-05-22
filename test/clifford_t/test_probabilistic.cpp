#include <gtest/gtest.h>

#include "type.hpp"
#include "su2.hpp"
#include "clifford_t/probabilistic.hpp"

using namespace::su2compiler;

// TEST(ProbSynthesize, ComputeOptProb) {
//     SU2 V = su2compiler::random_unitary();
//     int N = 50;
//     std::vector<SU2> availableU(N);
//     for(int i = 0; i < N; i++) availableU[i] = su2compiler::random_unitary();

//     for(int i = 0; i < N; i++){
//         std::cout << distance(V, availableU[i]) << std::endl;
//     }

//     auto [opt_dist, opt_prob] = clifford_t::compute_optimal_prob()
// }