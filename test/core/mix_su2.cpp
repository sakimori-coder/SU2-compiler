#include <gtest/gtest.h>

#include "core/type.hpp"
#include "core/su2.hpp"
#include "core/mix_su2.hpp"

using namespace su2compiler;

TEST(mix_su2Test, compute_optimal_prob)
{
    int N = 50;
    std::vector<SU2> availableU(N);
    for(int i = 0; i < N; i++) {
        availableU[i] = random_unitary();
    }

    MixSU2 mixing_unitary(availableU);

    SU2 targetV = random_unitary();
    math::sdp::Options opt;
    opt.MaxIteration = 200;
    opt.lambda       = 1e4;
    opt.betaStar     = 0.1;
    opt.betaBar      = 0.3;
    opt.gamma        = 0.9;
    opt.epsilon1     = 1e-30;
    opt.epsilon2     = 1e-30;
    opt.OUTPUT_HISTORY = true;
    auto res = mixing_unitary.compute_optimal_prob(targetV, opt);

    std::cout << res.obj_primal << std::endl;
    std::cout << res.obj_dual << std::endl;
    EXPECT_TRUE(res.status == STATUS::SUCCESS);

}