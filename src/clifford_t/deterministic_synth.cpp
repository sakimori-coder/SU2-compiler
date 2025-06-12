#include "clifford_t/deterministic_synth.hpp"

#include <string>
#include <vector>
#include <Eigen/Core>

#include "core/type.hpp"
#include "core/su2.hpp"
#include "math/constants.hpp"
#include "math/functions.hpp"
#include "math/lattice.hpp"
#include "math/linalg.hpp"
#include "clifford_t/exact_synth.hpp"
#include "util/time_profiler.hpp"

namespace su2compiler::clifford_t::deterministic {


using exact::U2Dzeta8;
template <typename T>
using Matrix8 = Eigen::Matrix<T, 8, 8>;
template <typename T>
using Vector8 = Eigen::Vector<T, 8>;


std::string synth(
        SU2 V,
        Real eps)
{
    int prec = std::max(53,
                        static_cast<int>(-8 * mpfr_get_exp(eps.mpfr_ptr()))
                       );
    Real::set_default_prec(prec);

    if(eps >= Real(1)) return "";

    Natural t = 0;
    while(true) {
        auto res = fixed_t_synth(V, eps, t);
        if(!res.empty()) {
            return exact::synth(res.front());
        }

        t++;
    }

    return "FAILURE";
}


std::vector<U2Dzeta8> fixed_t_synth(
        SU2 V,
        Real eps,
        Natural t)
{
    SU2 V_prime = V * SU2(math::ZETA16(), 
                          Complex(0));
    Natural t_prime = std::max(
                0,
                int(math::round(Real(int(t)) + 5.0/2.0 * math::log2(eps)).get_si())
    );

    std::vector<std::string> U_L_list;
    for(long mask = 0; mask < (1L << t_prime); mask++) {
        std::string sequence = "";
        for(Natural j = 0; j < t_prime; j++) {
            if(mask & (1L << j)) sequence += "HT";
            else                 sequence += "SHT";
        }
        U_L_list.push_back(sequence);
    }

    if(t_prime >= 1) {
        for(long mask = 0; mask < (1L << (t_prime - 1)); mask++) {
            std::string sequence = "T";
            for(Natural j = 0; j < t_prime-1; j++) {
                if(mask & (1L << j)) sequence += "HT";
                else                 sequence += "SHT";
            }
            U_L_list.push_back(sequence);
        }
    }

    std::vector<U2Dzeta8> U_list;
    for(auto U_L_seq : U_L_list) {
        U2Dzeta8 U_L(U_L_seq);

        if((t - t_prime) % 2 == 0) {
            Natural k = (t - t_prime + 2) / 2;
            auto U_R_list = solve_approx_lattice(
                    U_L.adjoint().toSU2() * V,
                    eps, 
                    k
            );
            for(auto U_R : U_R_list) {
                U_list.push_back(U_L * U2Dzeta8(U_R, 0, k));
            }
        } else {
            Natural k = (t - t_prime + 3) / 2;
            auto U_R_list = solve_approx_lattice(
                    U_L.adjoint().toSU2() * V_prime,
                    eps,
                    k
            );
            for(auto U_R : U_R_list) {
                U_list.push_back(U_L * U2Dzeta8(U_R, 1, k));
            }        
        }
    }

    std::erase_if(U_list, [&t](U2Dzeta8 U) {
        return exact::get_Tcount(U) != t; 
    });
    return U_list;
}


std::vector<ring::Zzeta8j> solve_approx_lattice(
        SU2 V,
        Real eps,
        Natural k)
{
    auto& prof = util::Profiler::instance();
    prof.start("SovleApproxLattice");

    Real r = math::pow_ui(math::SQRT2(), k);

    // std::cout << k << std::endl;
    Matrix8<Real> Sigma, Sigma_inv;
    Sigma << 1, math::INV_SQRT2(), 0,-math::INV_SQRT2(), 0, 0, 0, 0,
             0, math::INV_SQRT2(), 1, math::INV_SQRT2(), 0, 0, 0, 0,
             0, 0, 0, 0, 1, math::INV_SQRT2(), 0,-math::INV_SQRT2(),
             0, 0, 0, 0, 0, math::INV_SQRT2(), 1, math::INV_SQRT2(),
             1,-math::INV_SQRT2(), 0, math::INV_SQRT2(), 0, 0, 0, 0,
             0,-math::INV_SQRT2(), 1,-math::INV_SQRT2(), 0, 0, 0, 0,
             0, 0, 0, 0, 1,-math::INV_SQRT2(), 0, math::INV_SQRT2(),
             0, 0, 0, 0, 0,-math::INV_SQRT2(), 1,-math::INV_SQRT2();
    Sigma_inv = Sigma.transpose();
    Sigma_inv /= Real(2);
    
    Matrix8<Real> P;
    P << V.a,-V.b,-V.c,-V.d, 0, 0, 0, 0,
         V.b, V.a, V.d,-V.c, 0, 0, 0, 0,
         V.c,-V.d, V.a, V.b, 0, 0, 0, 0,
         V.d, V.c,-V.b, V.a, 0, 0, 0, 0,
         0, 0, 0, 0,         1, 0, 0, 0,
         0, 0, 0, 0,         0, 1, 0, 0,
         0, 0, 0, 0,         0, 0, 1, 0,
         0, 0, 0, 0,         0, 0, 0, 1; 

    Real e1 = r * (Real(1) - math::sqrt(Real(1) - eps*eps));
    e1 *= e1;
    Real e2 = r * eps;
    e2 *= e2;
    Real e3 = r;
    e3 *= e3;

    Matrix8<Real> Q = Matrix8<Real>::Zero();
    Q(0,0) = Real(1) / e1;
    Q(1,1) = Real(1) / e2;
    Q(2,2) = Real(1) / e2;
    Q(3,3) = Real(1) / e2;
    for(int i = 4; i < 8; i++) Q(i,i) = Real(1) / e3;
    Q /= Real(2);
    Q = Sigma.transpose() * P * Q * P.transpose() * Sigma;

    Vector8<Real> p;
    p << r * math::sqrt(Real(1) - eps*eps), 0, 0, 0, 0, 0, 0, 0;
    p = Sigma_inv * P * p;

    auto X = math::lattice::EnumIntegerPoints(Q, p, Real(1));

    Integer key = Integer(1) << k;
    std::vector<ring::Zzeta8j> Solutions;
    for(auto xvec : X) {
        ring::Zzeta8j q(xvec(0), xvec(1), xvec(2), xvec(3),
                        xvec(4), xvec(5), xvec(6), xvec(7));

        if(q.norm_quaternion() == key) {
            SU2 U(q.u.toComplex(), q.t.toComplex());
            U.unitalize();
            if(distance(U, V) < eps) {
                Solutions.push_back(q);
            }
        } 
    }

    prof.stop("SovleApproxLattice");
    return Solutions;
}




}