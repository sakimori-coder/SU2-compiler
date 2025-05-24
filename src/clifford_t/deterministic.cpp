#include "clifford_t/deterministic.hpp"

#include <algorithm>
#include <chrono>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/LU>

#include "ring/Zzeta8j.hpp"
#include "type.hpp"
#include "enum_integer_points.hpp"
#include "su2.hpp"
#include "clifford_t/exact_synthesis.hpp"


namespace su2compiler {
namespace clifford_t {


std::string deterministic_synthesis(SU2 V, Real eps) {
    if(eps >= Real(1, PrecBits)) return "";

    UINT t = 0;
    while(true) {
        std::cout << "t = " << t << std::endl;
        auto res = deterministic_synthesis_fixed_t(V, eps, t);
        if(!res.empty()) {
            std::cout << "START exact" << std::endl; 
            return exact_synthesis(res.front());
        }
            t++;
    }
    return "FAILURE";
}


std::vector<U2Dzeta8> deterministic_synthesis_fixed_t(SU2 V, Real eps, UINT t) {
    SU2 V_prime = V * SU2(ZETA16(), Complex(0.0, 0.0));
    UINT t_prime = std::max(0, int(round(Real(t, PrecBits) + 5.0/2.0*log2(eps))));
    
    std::vector<std::string> CT_sequence;
    for(UINT mask = 0; mask < (UINT(1)<<t_prime); mask++) {
        std::string sequence = "";
        for(int j = 0; j < t_prime; j++) {
            if(mask & (UINT(1)<<j)) sequence += "HT";
            else                    sequence += "SHT";
        }
        CT_sequence.push_back(sequence);
    }
    
    if(t_prime >= 1){
        for(UINT mask = 0; mask < (UINT(1)<<(t_prime-1)); mask++) {
            std::string sequence = "";
            for(int j = 0; j < t_prime-1; j++) {
                if(mask & (UINT(1)<<j)) sequence += "HT";
                else                    sequence += "SHT";
            }
            CT_sequence.push_back("T" + sequence);
        }
    }

    
    // std::cout << CT_sequence.size() << std::endl;
    std::vector<U2Dzeta8> U_list;
#pragma omp parallel for
    for(auto sequence : CT_sequence) {
        mpfr::mpreal::set_default_prec(PrecBits);
        
        U2Dzeta8 U_prime(sequence);

        int k_even = (t - t_prime + 3) / 2;
        int k_odd = (t - t_prime + 4) / 2;
        // std::cout << k_even << std::endl;
        // std::cout << k_odd << std::endl;

        std::chrono::system_clock::time_point start, end;
        double time;
        start = std::chrono::system_clock::now();
        if((t - t_prime) % 2 == 0) {
            UINT k_even = (t - t_prime + 3) / 2;
            auto U_even_list = solve_approx_synthesis(U_prime.adjoint().to_SU2() * V, eps, k_even);
#pragma omp critical
            for(auto U_even : U_even_list) U_list.push_back(U_prime * U2Dzeta8(U_even, 0, k_even));
        } else {
            UINT k_odd = (t - t_prime + 4) / 2;
            auto U_odd_list = solve_approx_synthesis(U_prime.adjoint().to_SU2() * V_prime, eps, k_odd);
#pragma omp critical
            for(auto U_odd : U_odd_list) U_list.push_back(U_prime * U2Dzeta8(U_odd, 1, k_odd));
        }
        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        // std::cout << time << "[ms]" << std::endl;
    }

    std::erase_if(U_list, [&t](U2Dzeta8 U) { return get_Tcount(U) != t; });
    return U_list;
}


std::vector<ring::Zzeta8j> solve_approx_synthesis(SU2 V, Real eps, UINT k) {
    Real r = pow_ui(SQRT2(), k);

    // std::cout << k << std::endl;
    Matrix8R Sigma, Sigma_inv;
    Sigma << 1, INV_SQRT2(), 0,-INV_SQRT2(), 0, 0, 0, 0,
             0, INV_SQRT2(), 1, INV_SQRT2(), 0, 0, 0, 0,
             0, 0, 0, 0, 1, INV_SQRT2(), 0,-INV_SQRT2(),
             0, 0, 0, 0, 0, INV_SQRT2(), 1, INV_SQRT2(),
             1,-INV_SQRT2(), 0, INV_SQRT2(), 0, 0, 0, 0,
             0,-INV_SQRT2(), 1,-INV_SQRT2(), 0, 0, 0, 0,
             0, 0, 0, 0, 1,-INV_SQRT2(), 0, INV_SQRT2(),
             0, 0, 0, 0, 0,-INV_SQRT2(), 1,-INV_SQRT2();
    Sigma_inv = Sigma.transpose();
    Sigma_inv /= 2.0;
    
    Matrix8R D = Matrix8R::Zero();
    Matrix8R P;
    P << V.a,-V.b,-V.c,-V.d, 0, 0, 0, 0,
         V.b, V.a, V.d,-V.c, 0, 0, 0, 0,
         V.c,-V.d, V.a, V.b, 0, 0, 0, 0,
         V.d, V.c,-V.b, V.a, 0, 0, 0, 0,
         0, 0, 0, 0,         1, 0, 0, 0,
         0, 0, 0, 0,         0, 1, 0, 0,
         0, 0, 0, 0,         0, 0, 1, 0,
         0, 0, 0, 0,         0, 0, 0, 1; 
    Real e1 = r * (Real(1, PrecBits) - sqrt(Real(1, PrecBits) - eps*eps));
    e1 *= e1;
    Real e2 = r * eps;
    e2 *= e2;
    Real e3 = r;
    e3 *= e3;

    D(0,0) = e1;
    D(1,1) = e2;
    D(2,2) = e2;
    D(3,3) = e2;
    for(int i = 4; i < 8; i++) D(i,i) = e3;
    D *= Real(2, PrecBits);
    
    D = Sigma_inv * P * D * P.transpose() * Sigma_inv.transpose();
    Matrix8R Q = D.inverse();

    Vector8R p;
    p << r * sqrt(1 - eps*eps), 0, 0, 0, 0, 0, 0, 0;
    p = Sigma_inv * P * p;

    std::cout << "START ENUM" << std::endl;
    auto X = EnumIntegerPoints(Q, p, Real(1, PrecBits));
    std::cout << "END ENUM" << std::endl;
    
    Integer key = Integer(1) << k;
    std::vector<ring::Zzeta8j> Solutions;
    for(auto xvec : X) {
        ring::Zzeta8j q(xvec(0), xvec(1), xvec(2), xvec(3),
                        xvec(4), xvec(5), xvec(6), xvec(7));

        if(q.norm_quaternion() == key) {
            SU2 U(q.u.to_Complex(), q.t.to_Complex());
            U.unitalize();
            if(distance(U, V) < eps) {
                Solutions.push_back(q);
            }
        } 
    }

    // auto divisibleBySqrt2 = [](ring::Zzeta8j q) { return q.divisible(ring::Zroot2(0,1)); };
    // Solutions.erase(std::remove_if(Solutions.begin(), Solutions.end(), divisibleBySqrt2), Solutions.end());

    return Solutions;
}


}
}