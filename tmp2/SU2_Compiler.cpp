#include "su2compiler.hpp"

#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <Eigen/Core>
#include <Eigen/LU>
#include <chrono>

#include "type.hpp"
#include "rings.hpp"
#include "SpecialUnitary2.hpp"
#include "Clifford_T_1Q.hpp"
#include "enum_integer_points.hpp"



namespace su2compiler
{


std::string Deterministic_Clifford_T(SpecialUnitary2 V, Real eps)
{
    std::vector<Clifford_T_1Q> candidates;

    // case T-count is even
    auto solutions_even = solve_approximation_synthesis_CliffordT(V, eps);
    for(auto [u, t, k] : solutions_even){
        candidates.push_back(Clifford_T_1Q(u, t, 0, k));
    }

    // case T-count is odd
    SpecialUnitary2 V_prime;
    V_prime = SpecialUnitary2(sqrt_omega.real(), sqrt_omega.imag(), 0, 0) * V;
    auto solutions_odd = solve_approximation_synthesis_CliffordT(V_prime, eps);
    for(auto [u, t, k] : solutions_odd){
        candidates.push_back(Clifford_T_1Q(u, t, 1, k));
    }

    Clifford_T_1Q U = *std::min_element(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
        return a.compute_Tcount() < b.compute_Tcount();
    });
    
    return U.compute_sequence();
}  



// std::vector<std::tuple<ZOmega, ZOmega, Integer>> solve_approximation_synthesis_CliffordT(SpecialUnitary2 V, Real eps)
// {
//     Matrix8R P;
//     P << V.a,-V.b,-V.c,-V.d, 0, 0, 0, 0,
//          V.b, V.a, V.d,-V.c, 0, 0, 0, 0,
//          V.c,-V.d, V.a, V.b, 0, 0, 0, 0,
//          V.d, V.c,-V.b, V.a, 0, 0, 0, 0,
//          0, 0, 0, 0,         1, 0, 0, 0,
//          0, 0, 0, 0,         0, 1, 0, 0,
//          0, 0, 0, 0,         0, 0, 1, 0,
//          0, 0, 0, 0,         0, 0, 0, 1;  

//     Matrix8R Sigma, Sigma_inv;
//     Sigma << 1, 1/sqrt2, 0,-1/sqrt2, 0, 0, 0, 0,
//              0, 1/sqrt2, 1, 1/sqrt2, 0, 0, 0, 0,
//              0, 0, 0, 0, 1, 1/sqrt2, 0,-1/sqrt2,
//              0, 0, 0, 0, 0, 1/sqrt2, 1, 1/sqrt2,
//              1,-1/sqrt2, 0, 1/sqrt2, 0, 0, 0, 0,
//              0,-1/sqrt2, 1,-1/sqrt2, 0, 0, 0, 0,
//              0, 0, 0, 0, 1,-1/sqrt2, 0, 1/sqrt2,
//              0, 0, 0, 0, 0,-1/sqrt2, 1,-1/sqrt2;

//     Sigma_inv = Sigma.transpose();
//     Sigma_inv /= 2.0;

//     int k = 0;
//     while(true){
//         std::cout << "k=" << k << std::endl;
//         Real r = pow(sqrt2, (Real)k);
//         std::vector<std::tuple<ZOmega, ZOmega, Integer>> solutions;

//         Matrix8R D = Matrix8R::Zero();
//         Real e1 = r * (1 - sqrt(1 - eps*eps));
//         e1 *= e1;
//         Real e2 = r * eps;
//         e2 *= e2;
//         Real e3 = r;
//         e3 *= e3;

//         D(0,0) = e1;
//         D(1,1) = e2;
//         D(2,2) = e2;
//         D(3,3) = e2;
//         for(int i = 4; i < 8; i++) D(i,i) = e3;
//         D *= (Real)2.0;
        
//         D = Sigma_inv * P * D * P.transpose() * Sigma_inv.transpose();
//         Matrix8R Q = D.inverse();

//         Vector8R p;
//         p << r * sqrt(1 - eps*eps), 0, 0, 0, 0, 0, 0, 0;
//         p = Sigma_inv * P * p;

//         std::chrono::system_clock::time_point start, end;
//         double time;
//         start = std::chrono::system_clock::now();
//         auto Xvec = EnumIntegerPoints(Q, p, 1.0);
//         std::cout << "|X|=" << Xvec.size() << std::endl;
//         end = std::chrono::system_clock::now();
//         time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
//         std::cout << time << "[ms]" << std::endl;

//         start = std::chrono::system_clock::now();
//         Integer key = (Integer)1 << k;
//         for(auto xvec : Xvec){
//             ZOmega u(xvec(0), xvec(1), xvec(2), xvec(3));
//             ZOmega t(xvec(4), xvec(5), xvec(6), xvec(7));

//             if(u * complex_conj(u) + t * complex_conj(t) == key){
//                 SpecialUnitary2 U(u.to_Complex().real(),
//                                   u.to_Complex().imag(),
//                                   t.to_Complex().real(),
//                                   t.to_Complex().imag());
//                 U.unitalize();
//                 if(distance(U, V) < eps){
//                     solutions.push_back({u, t, k});
//                 }
//             }
//         }
//         end = std::chrono::system_clock::now();
//         time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
//         std::cout << time << "[ms]" << std::endl;

//         if(!solutions.empty()) return solutions;
//         k++;
//     }
// }



std::vector<std::tuple<ZOmega, ZOmega, Integer>> solve_approximation_synthesis_CliffordT(SpecialUnitary2 V, Real eps)
{
    Matrix8R Sigma, Sigma_inv;
    Sigma << 1, 1/sqrt2, 0,-1/sqrt2, 0, 0, 0, 0,
             0, 1/sqrt2, 1, 1/sqrt2, 0, 0, 0, 0,
             0, 0, 0, 0, 1, 1/sqrt2, 0,-1/sqrt2,
             0, 0, 0, 0, 0, 1/sqrt2, 1, 1/sqrt2,
             1,-1/sqrt2, 0, 1/sqrt2, 0, 0, 0, 0,
             0,-1/sqrt2, 1,-1/sqrt2, 0, 0, 0, 0,
             0, 0, 0, 0, 1,-1/sqrt2, 0, 1/sqrt2,
             0, 0, 0, 0, 0,-1/sqrt2, 1,-1/sqrt2;

    Sigma_inv = Sigma.transpose();
    Sigma_inv /= 2.0;

    int k = 0;
    while(true){
        std::cout << "k=" << k << std::endl;
        Real r = pow(sqrt2, (Real)k);
        std::vector<std::tuple<ZOmega, ZOmega, Integer>> solutions;

        Real delta = min((Real)1.0, 2.0*pow(r, -(Real)1.6) / eps);
        int N = (int)ceil(1.0/ sqrt(8.0/(PI*PI) * (1 - pow(1-delta*delta, 1.0/(Real)3.0))));

        std::vector<Real> theta1_list(N), theta2_list(N), theta3_list(2*N);
        for(int i = 0; i < N; i++) theta1_list[i] = PI * ((Real)i / N);
        for(int i = 0; i < N; i++) theta2_list[i] = PI * ((Real)i / N);
        for(int i = 0; i < 2*N; i++) theta3_list[i] = PI * ((Real)i / N);

#pragma omp parallel for collapse(3)
        for(auto &theta1 : theta1_list){
            for(auto &theta2 : theta2_list){
                for(auto &theta3 : theta3_list){
                    std::chrono::system_clock::time_point start, end;
                    double time;
                    start = std::chrono::system_clock::now();

                    Real qx = cos(theta1);
                    Real qy = sin(theta1) * cos(theta2);
                    Real qz = sin(theta1) * sin(theta2) * cos(theta3);
                    Real qw = sin(theta1) * sin(theta2) * sin(theta3);

                    Matrix8R P;
                    P << V.a,-V.b,-V.c,-V.d, 0, 0, 0, 0,
                         V.b, V.a, V.d,-V.c, 0, 0, 0, 0,
                         V.c,-V.d, V.a, V.b, 0, 0, 0, 0,
                         V.d, V.c,-V.b, V.a, 0, 0, 0, 0,
                         0, 0, 0, 0, qx,-qy,-qz,-qw,
                         0, 0, 0, 0, qy, qx, qw,-qz,
                         0, 0, 0, 0, qz,-qw, qx, qy,
                         0, 0, 0, 0, qw, qz,-qy, qx;

                    Matrix8R D = Matrix8R::Zero();
                    Real e1 = r * (1 - sqrt(1 - eps*eps));
                    e1 *= e1;
                    Real e2 = r * eps;
                    e2 *= e2;
                    Real e3 = r * (delta*delta);
                    e3 *= e3;
                    Real e4 = r * sqrt(1 - (1-delta*delta)*(1-delta*delta));
                    e4 *= e4;

                    D(0,0) = e1;
                    for(int i = 1; i < 4; i++) D(i,i) = e2;
                    D(4,4) = e3;
                    for(int i = 5; i < 8; i++) D(i,i) = e4;
                    D *= (Real)2.0;

                    D = Sigma_inv * P * D * P.transpose() * Sigma_inv.transpose();
                    Matrix8R Q = D.inverse();

                    Vector8R p;
                    p << r * sqrt(1 - eps*eps), 0, 0, 0,
                         r * (1 - delta*delta), 0, 0, 0;
                    p = Sigma_inv * P * p;

                    auto Xvec = EnumIntegerPoints(Q, p, 1.0);
                    Integer key = (Integer)1 << k;
                    for(auto xvec : Xvec){
                        ZOmega u(xvec(0), xvec(1), xvec(2), xvec(3));
                        ZOmega t(xvec(4), xvec(5), xvec(6), xvec(7));

                        if(u * complex_conj(u) + t * complex_conj(t) == key){
                            SpecialUnitary2 U(u.to_Complex().real(),
                                            u.to_Complex().imag(),
                                            t.to_Complex().real(),
                                            t.to_Complex().imag());
                            U.unitalize();
                            if(distance(U, V) < eps){
#pragma omp critical
                                solutions.push_back({u, t, k});
                            }
                        }
                    }
                    end = std::chrono::system_clock::now();
                    time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
                    std::cout << "Total:" << time << "[ms]" << std::endl;

                }
            }
        }
        std::cout << "分割数 : " << theta1_list.size() * theta2_list.size() * theta3_list.size() << std::endl;
        if(!solutions.empty()) return solutions;
        k++;
    }
}


}