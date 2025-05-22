#include <gtest/gtest.h>

#include <complex>
#include <string>

#include "type.hpp"
#include "ring/all.hpp"
#include "su2.hpp"
#include "clifford_t/deterministic.hpp"

using namespace::su2compiler;

TEST(deterministic_synthesis, random) {
    SU2 V = su2compiler::random_unitary(1234);
    // SU2 V(Real(3./5.), 0.0, Real(4./5.), 0.0);
    Real eps = 1e-5;

    std::cout << "V = \n" << V << std::endl;

    std::string sequence = clifford_t::deterministic_synthesis(V, eps);

    SU2 H(Complex(0.0, 1.0) / SQRT2(), Complex(0.0, 1.0) / SQRT2());
    SU2 S(std::conj(ZETA8()), Complex(0.0, 0.0));
    SU2 T(std::conj(ZETA16()), Complex(0.0, 0.0));

    SU2 U(1,0,0,0);
    for(char gate : sequence) {
        switch (gate)
        {
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

    double dist = static_cast<double>(distance(U,V));
    EXPECT_NEAR(dist, 0.0, static_cast<double>(eps));
    std::cout << V << std::endl;
    std::cout << sequence << std::endl;
    // std::cout << "#T=" << std::count(sequence.begin(), sequence.end(), 'T') << std::endl;
}


// TEST(SolveApproxSynthesis, example) {
//     SU2 V = random_unitary(1234);
//     std::cout << V << std::endl;
//     Real eps = 1e-5;
//     int k = 25;
//     auto U = clifford_t::solve_approx_synthesis(V, eps, k);
//     std::set<ring::Zzeta8j> Us(U.begin(), U.end());
//     for(auto q : Us) {
//         // std::cout << q << std::endl;
//         std::cout << (q.u.a & 1) << (q.u.b & 1) << (q.u.c & 1) << (q.u.d & 1) << std::endl;
//         std::cout << (q.t.a & 1) << (q.t.b & 1) << (q.t.c & 1) << (q.t.d & 1) << std::endl;
//         std::cout << std::endl;
//     }

//     ring::Zzeta8j S(0,1,0,0, 0,0,0,0);
//     ring::Zzeta8j H(1,0,0,0, 1,0,0,0);
//     ring::Zzeta8j Hz(1,0,0,0, 0,1,0,0);

//     static const ring::Zzeta8j W[6] = {
//         S * Hz * S * H,
//         S * Hz * H,
//         Hz * S * H,
//         Hz * H,
//         S * H * ring::Zroot2(0,1),
//         H * ring::Zroot2(0,1)
//     };

    
//     for(auto q : Us) {
//         // std::cout << (q.u.a & 1) << (q.u.b & 1) << (q.u.c & 1) << (q.u.d & 1) << std::endl;
//         // std::cout << (q.t.a & 1) << (q.t.b & 1) << (q.t.c & 1) << (q.t.d & 1) << std::endl;
//         std::cout << q << std::endl;
//         std::cout << q.mod2() << std::endl;
//         q *= ring::Zroot2(0,1);
//         for(auto w : W) {
//             if(q.leftDivisible(w)) std::cout << "YES" << std::endl;
//             else std::cout << "NO" << std::endl;
//         }
//         // if(q.divisible(ring::Zroot2(0,1))) std::cout << "DIVISIBLE √2" << std::endl;
//         std::cout << std::endl;
//     }

//     for(auto w : W){
//         std::cout << w.norm_quaternion() << std::endl;
//     }
// }



// TEST(tmp, tmp) {
//     using namespace std;
//     using namespace ring;
//     Zzeta8j q(2391, 610, 645, -1563, 4703, -348, 1350, 738);
//     q *= ring::Zzeta8j(0,0,1,0, 0,0,0,0);

//     ring::Zzeta8j S(0,1,0,0, 0,0,0,0);
//     ring::Zzeta8j H(1,0,0,0, 1,0,0,0);
//     ring::Zzeta8j Hz(1,0,0,0, 0,1,0,0);

//     q  = S * q;

//     cout << q << endl;
//     cout << q.mod2() << endl;

//     q = Hz * q;
//     cout << q << endl;
//     q /= Zroot2(0,1);
//     cout << q << std::endl;
//     cout << q.mod2() << endl;

//     cout << q.norm_quaternion() << endl;
// }
