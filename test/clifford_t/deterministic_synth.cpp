#include <gtest/gtest.h>

#include "core/type.hpp"
#include "core/su2.hpp"
#include "ring/all.hpp"
#include "clifford_t/deterministic_synth.hpp"

using namespace su2compiler;



TEST(CliffordTDetSynthTest, random) {
    SU2 V = random_unitary(1234);
    Real eps = 1e-5;

    std::string sequence = clifford_t::deterministic::synth(V, eps);
    // std::cout << "T-count = " << std::count(sequence.begin(), sequence.end(), 'T') << std::endl;

    SU2 H(Complex(0, 1) / math::SQRT2(), 
          Complex(0, 1) / math::SQRT2());
    SU2 S(std::conj(math::ZETA8()), Complex(0, 0));
    SU2 T(std::conj(math::ZETA16()), Complex(0, 0));

    SU2 U(1,0,0,0);
    for(char gate : sequence) {
        switch (gate) {
            case 'H':
                U *= H;
                break;
            case 'S':
                U *= S;
                break;
            case 'T':
                U *= T;
                break;
            default:
                break;
        }
    }

    Real dist = distance(U, V);
    EXPECT_TRUE(dist < eps);
}