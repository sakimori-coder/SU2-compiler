#include <gtest/gtest.h>

#include <complex>
#include <string>

#include "type.hpp"
#include "ring/all.hpp"
#include "su2.hpp"
#include "clifford_t/synthesize.hpp"

using namespace::su2_compiler;

// TEST(Snthesize, random) {
//     SU2 V = su2_compiler::random_unitary(1234);
//     Real eps = 1e-7;

//     std::string sequence = clifford_t::synthesize(V, eps);
    
//     SU2 H(Complex(0.0, 1.0) / SQRT2, Complex(0.0, 1.0) / SQRT2);
//     SU2 S(std::conj(ZETA8), Complex(0.0, 0.0));
//     SU2 T(std::conj(ZETA16), Complex(0.0, 0.0));

//     SU2 U(1,0,0,0);
//     for(char gate : sequence) {
//         switch (gate)
//         {
//         case 'H':
//             U *= H;
//             break;
//         case 'S':
//             U *= S;
//             break;
//         case 'T':
//             U *= T;
//             break;
//         default:
//             break;
//         }
//     }

//     double dist = static_cast<double>(distance(U,V));
//     EXPECT_NEAR(dist, 0.0, static_cast<double>(eps));
// }


TEST(SolveApproxSynthesis, example) {
    SU2 V = random_unitary(1234);
    std::cout << V << std::endl;
    Real eps = 1e-5;
    int k = 25;
    auto U = clifford_t::solve_approx_synthesis(V, eps, k);
    std::set<ring::Zzeta8j> Us(U.begin(), U.end());
    for(auto q : Us) {
        std::cout << q << std::endl;
    }
    static const ring::Zzeta8j W[4] = {
        ring::Zzeta8j(1,0,0,0, 1,0,0,0),
        ring::Zzeta8j(1,0,0,0, 0,1,0,0),
        ring::Zzeta8j(1,0,0,0, 0,0,1,0),
        ring::Zzeta8j(1,0,0,0, 0,0,0,1),
    };

    
    for(auto q : Us) {
        std::cout << q << std::endl;
        for(auto w : W) {
            if(q.leftDivisible(w)) std::cout << "YES" << std::endl;
            else std::cout << "NO" << std::endl;
        }
        if(q.divisible(ring::Zroot2(0,1))) std::cout << "DIVISIBLE √2" << std::endl;
        std::cout << std::endl;
    }
}