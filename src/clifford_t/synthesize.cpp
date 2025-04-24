#include "clifford_t/synthesize.hpp"

#include <algorithm>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/LU>

#include "ring/Zzeta8j.hpp"
#include "type.hpp"
#include "enum_integer_points.hpp"
#include "su2.hpp"
#include "clifford_t/exact_synthesize.hpp"


namespace su2_compiler {
namespace clifford_t {

std::string synthesize(SU2 V, Real eps) {
    SU2 V_prime = V * SU2(ZETA16, Complex(0.0, 0.0));

    int k = 0;
    while(true) {
        std::cout << "k=" << k << std::endl;
        // std::cout << "start" << std::endl;
        auto U_even_list = solve_approx_synthesis(V, eps, k);
        // std::cout << "end" << std::endl;
        auto U_odd_list = solve_approx_synthesis(V_prime, eps, k);

        std::vector<U2Dzeta8> U_list;
        for(auto U_even : U_even_list) U_list.push_back(U2Dzeta8(U_even, 0, k));
        // for(auto U_odd : U_odd_list) U_list.push_back(U2Dzeta8(U_odd, 1, k));

        if(!U_list.empty()) {
            U2Dzeta8 U = *std::min_element(U_list.begin(), U_list.end(), [](auto a, auto b) {
                return get_Tcount(a) < get_Tcount(b);
            });
            return exact_synthesize(U);
        }
        k++;
    }
}


static const ring::Zzeta8j W[4] = {
    ring::Zzeta8j(1,0,0,0, 1,0,0,0),
    ring::Zzeta8j(1,0,0,0, 0,1,0,0),
    ring::Zzeta8j(1,0,0,0, 0,0,1,0),
    ring::Zzeta8j(1,0,0,0, 0,0,0,1),
};

static const ring::Zzeta8j W_conj[4] = {
    W[0].conj_quaternion(),
    W[1].conj_quaternion(),
    W[2].conj_quaternion(),
    W[3].conj_quaternion(),
};

std::vector<ring::Zzeta8j> solve_approx_synthesis(SU2 V, Real eps, int k) {

    Real r = pow(SQRT2, Real(k));
    
    const Real threshold = 1e3;   // hyper parameter
    // if(r*r*r*r * r*r*r*r * eps*eps*eps*eps*eps >= threshold){
    if(k == 26) {
        std::vector<ring::Zzeta8j> Solutions;
// #pragma omp parallel for
        for(int i = 0; i < 4; i++) {
            std::cout << SU2(W_conj[i].u.to_Complex() / SQRT2, W_conj[i].t.to_Complex() / SQRT2) << std::endl; 
            auto Solutions_sub = solve_approx_synthesis(
                    SU2(W_conj[i].u.to_Complex() / SQRT2, W_conj[i].t.to_Complex() / SQRT2) * V, eps, k-1);
// #pragma omp critical
            for(auto U : Solutions_sub) Solutions.push_back(W[i] * U);
            for(auto U : Solutions_sub) {
                std::cout << W[i] * U << std::endl;
            }
            std::cout << std::endl;
        }
        return Solutions;
    }

    // std::cout << k << std::endl;
    Matrix8R Sigma, Sigma_inv;
    Sigma << 1, 1/SQRT2, 0,-1/SQRT2, 0, 0, 0, 0,
             0, 1/SQRT2, 1, 1/SQRT2, 0, 0, 0, 0,
             0, 0, 0, 0, 1, 1/SQRT2, 0,-1/SQRT2,
             0, 0, 0, 0, 0, 1/SQRT2, 1, 1/SQRT2,
             1,-1/SQRT2, 0, 1/SQRT2, 0, 0, 0, 0,
             0,-1/SQRT2, 1,-1/SQRT2, 0, 0, 0, 0,
             0, 0, 0, 0, 1,-1/SQRT2, 0, 1/SQRT2,
             0, 0, 0, 0, 0,-1/SQRT2, 1,-1/SQRT2;
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
    Real e1 = r * (1 - sqrt(1 - eps*eps));
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
    D *= (Real)2.0;
    
    D = Sigma_inv * P * D * P.transpose() * Sigma_inv.transpose();
    Matrix8R Q = D.inverse();

    Vector8R p;
    p << r * sqrt(1 - eps*eps), 0, 0, 0, 0, 0, 0, 0;
    p = Sigma_inv * P * p;

    auto X = EnumIntegerPoints(Q, p, 1.0);
    
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

    auto divisibleBySqrt2 = [](ring::Zzeta8j q) { return q.divisible(ring::Zroot2(0,1)); };
    // Solutions.erase(std::remove_if(Solutions.begin(), Solutions.end(), divisibleBySqrt2), Solutions.end());

    return Solutions;
}


}
}