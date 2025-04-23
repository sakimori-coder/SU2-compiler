#include "Prob_Synthesis.hpp"

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <Eigen/Core>

#include "type.hpp"
#include "rings.hpp"
#include "quaternion.hpp"
#include "U2_ZOmega.hpp"
#include "ExactSynthesis.hpp"
#include "eps_net_verification.hpp"
#include "sdpa_dd.hpp"
#include "SU2_compiler.cpp"


namespace SU2_Compiler
{
    Eigen::Matrix<FTYPE, 4, 4> CJ_MB(Eigen::Matrix<CTYPE, 2, 2> U)
    {
        Eigen::Matrix<CTYPE, 4, 4> CJ_U = Choi_Jamiolkowski<2>(U);

        CTYPE i(0., 1.);
        CTYPE iinv_sqrt2(0, inv_sqrt2);
        Eigen::Matrix<CTYPE, 4, 4> M = Eigen::Matrix<CTYPE, 4, 4>::Zero(4,4);
        M(0,0) =  inv_sqrt2;
        M(0,3) =  inv_sqrt2;
        M(1,0) =-iinv_sqrt2;
        M(1,3) = iinv_sqrt2;
        M(2,1) =-iinv_sqrt2;
        M(2,2) =-iinv_sqrt2;
        M(3,1) =  inv_sqrt2;
        M(3,2) =- inv_sqrt2;

        Eigen::Matrix<CTYPE, 4, 4> CJ_U_MB;
        CJ_U_MB = M * CJ_U * M.adjoint();

        return CJ_U_MB.real();
    }


    template <typename T, int N, int M>
    std::vector<std::vector<T>> to_2d_vec(const Eigen::Matrix<T, N, M>& A)
    {
        std::vector<std::vector<T>> ret(N, std::vector<T>(M));
        
        for(int i = 0; i < N; i++){
            for(int j = 0; j < M; j++) ret[i][j] = A(i,j);
        }
        return ret;
    }


    FTYPE distance(quaternion A, std::vector<quaternion> B, std::vector<FTYPE> prob)
    {
        const int N = B.size();
        // std::cout << "N " << N << std::endl;
        
        Eigen::Matrix<FTYPE, 4, 4> CJ_A_MB = CJ_MB(A.get_Matrix());


        std::vector<Eigen::Matrix<FTYPE, 4, 4>> CJ_B_MB(N);
        for(int i = 0; i < N; i++) CJ_B_MB[i] = CJ_MB(B[i].get_Matrix());

        std::vector<Eigen::Matrix<FTYPE, 4, 4>> Symmetric_basis(10, Eigen::Matrix<FTYPE, 4, 4>::Zero(4,4));
        for(int i = 0; i < 4; i++) Symmetric_basis[i](i,i) = 1;
        Symmetric_basis[4](0,1) = Symmetric_basis[4](1,0) = 1;
        Symmetric_basis[5](0,2) = Symmetric_basis[5](2,0) = 1;
        Symmetric_basis[6](0,3) = Symmetric_basis[6](3,0) = 1;
        Symmetric_basis[7](1,2) = Symmetric_basis[7](2,1) = 1;
        Symmetric_basis[8](1,3) = Symmetric_basis[8](3,1) = 1;
        Symmetric_basis[9](2,3) = Symmetric_basis[9](3,2) = 1;

        int num_variable = 10;

        std::vector<FTYPE> C(num_variable, 0);
        C[0] = C[1] = C[2] = C[3] = 0.5;

        std::vector<std::vector<FTYPE>> ZERO_mat(4, std::vector<FTYPE>(4, 0));
        vector_4d<FTYPE> F(num_variable + 1, std::vector<std::vector<std::vector<FTYPE>>>(2));
        
        F[0][0] = ZERO_mat;
        Eigen::Matrix<FTYPE, 4, 4> D = -CJ_A_MB;
        for(int i = 0; i < N; i++) D += prob[i]*CJ_B_MB[i];
        F[0][1] = to_2d_vec<FTYPE, 4, 4>(D);

        for(int i = 1; i <= 10; i++){
            F[i][0] = to_2d_vec<FTYPE, 4, 4>(Symmetric_basis[i-1]);
            F[i][1] = to_2d_vec<FTYPE, 4, 4>(Symmetric_basis[i-1]);
        } 

        std::vector<FTYPE> xVec = sdpa_solver(C, F);

        return (xVec[0] + xVec[1] + xVec[2] + xVec[3]) / 2;
    }

    FTYPE distance(quaternion U, std::vector<std::pair<FTYPE, U2_ZOmega>> Prob_Clifford_T)
    {
        std::vector<quaternion> B;
        std::vector<FTYPE> prob;
        for(auto [p, U_ZOmega] : Prob_Clifford_T){
            B.push_back(to_quaternion(U_ZOmega));
            prob.push_back(p);
        }    
        return distance(U, B, prob);
    }

    FTYPE distance(quaternion U, std::vector<std::pair<FTYPE, std::string>> Prob_Clifford_T)
    {
        std::vector<quaternion> B;
        std::vector<FTYPE> prob;
        for(auto [p, seq] : Prob_Clifford_T){
            B.push_back(to_quaternion(seq));
            prob.push_back(p);
        }
        return distance(U, B, prob);
    }

    std::pair<FTYPE, std::vector<FTYPE>> get_optimal_prob(std::vector<quaternion> availableU, quaternion targetU)
    {
        const int N = availableU.size();
        
        std::vector<Eigen::Matrix<FTYPE, 4, 4>> availableCJUMB(N);
        for(int i = 0; i < N; i++) availableCJUMB[i] = CJ_MB(availableU[i].get_Matrix());
        Eigen::Matrix<FTYPE, 4, 4> targetCJUMB = CJ_MB(targetU.get_Matrix());
        
        
        std::vector<Eigen::Matrix<FTYPE, 4, 4>> Symmetric_basis(10, Eigen::Matrix<FTYPE, 4, 4>::Zero(4,4));
        for(int i = 0; i < 4; i++) Symmetric_basis[i](i,i) = 1;
        Symmetric_basis[4](0,1) = Symmetric_basis[4](1,0) = 1;
        Symmetric_basis[5](0,2) = Symmetric_basis[5](2,0) = 1;
        Symmetric_basis[6](0,3) = Symmetric_basis[6](3,0) = 1;
        Symmetric_basis[7](1,2) = Symmetric_basis[7](2,1) = 1;
        Symmetric_basis[8](1,3) = Symmetric_basis[8](3,1) = 1;
        Symmetric_basis[9](2,3) = Symmetric_basis[9](3,2) = 1;
        
        int num_variable = 10 + N;

        // 変数の並び順は (Sの変数10個, 混合確率)
        std::vector<FTYPE> C(num_variable, 0);
        C[0] = C[1] = C[2] = C[3] = 1 / (FTYPE)0.5;

        std::vector<std::vector<FTYPE>> ZERO_mat(4, std::vector<FTYPE>(4, 0));
        std::vector<std::vector<FTYPE>> ZERO_vec(1, std::vector<FTYPE>(N, 0)); 
        vector_4d<FTYPE> F(num_variable + 1, std::vector<std::vector<std::vector<FTYPE>>>(4));
        
        // F0
        F[0][0] = ZERO_mat;
        F[0][1] = to_2d_vec<FTYPE, 4, 4>(targetCJUMB);
        F[0][2] = ZERO_vec;
        F[0][3] = {{-1}};

        // F1~F10
        for(int i = 1; i <= 10; i++){
            F[i][0] = to_2d_vec<FTYPE, 4, 4>(Symmetric_basis[i-1]);
            F[i][1] = to_2d_vec<FTYPE, 4, 4>(Symmetric_basis[i-1]);
            F[i][2] = ZERO_vec;
            F[i][3] = {{0}};
        }

        // F11~FN
        for(int i = 0; i < N; i++){
            F[i+11][0] = ZERO_mat;
            F[i+11][1] = to_2d_vec<FTYPE, 4, 4>(availableCJUMB[i]);
            F[i+11][2] = ZERO_vec; 
            F[i+11][2][0][i] = 1;
            F[i+11][3] = {{-1}};
        }
        
        std::vector<FTYPE> xVec = sdpa_solver(C, F);    
        std::vector<FTYPE> prob(N);
        for(int i = 0; i < N; i++) prob[i] = xVec[i+10];

        FTYPE opt_val = (xVec[0] + xVec[1] + xVec[2] + xVec[3]) / 2;

        return {opt_val, prob};
    }


    // 与えられたTカウント以下で作れる最適な確率混合ユニタリを求める。もし、下限が目的のepsよりも大きい場合は計算を行わない。(本当は奇数・偶数で前の結果を再利用できるが、可読性が著しく悪いので行わない)
    std::vector< std::pair<FTYPE, U2_ZOmega>>
    optimal_prob_unitary(int T_count, quaternion targetU, FTYPE target_eps)
    {
        int k_odd, k_even;
        if(T_count % 2){
            k_odd = (T_count + 3) / 2;
            k_even = (T_count + 1) / 2;
        }else{
            k_odd = (T_count + 2) / 2;
            k_even = (T_count + 2) / 2;
        }

        FTYPE eps_prime = pow(2.0, -(FTYPE)T_count / 3.0) / 2.0;   // T_count = 3log_2(1/2eps')をeps'について解いた
        // FTYPE eps_prime = sqrt(target_eps) / 2;
        FTYPE delta_eps_prime = eps_prime / 5;


        while(true){
            std::vector<U2_ZOmega> availableU_odd = enum_u_t(targetU, 2.0*eps_prime, k_odd, 1);
            std::vector<U2_ZOmega> availableU_even = enum_u_t(targetU, 2.0*eps_prime, k_even, 0);
            std::vector<quaternion> availableU;
            std::vector<U2_ZOmega> availableU_ZOmega;

            for(auto &U_ZOmega : availableU_odd){
                if(get_T_count(U_ZOmega) > T_count) continue;

                quaternion U = to_quaternion(U_ZOmega);
                availableU.push_back(U);
                availableU_ZOmega.push_back(U_ZOmega);
                //std::cout << U << std::endl;
            }

            for(auto &U_ZOmega : availableU_even){
                if(get_T_count(U_ZOmega) > T_count) continue;

                quaternion U = to_quaternion(U_ZOmega);
                availableU.push_back(U);
                availableU_ZOmega.push_back(U_ZOmega);
                // std::cout << U << std::endl;
            }

            FTYPE lower = 1.0;   // availableUで作ることができる確率混合ユニタリの誤差の下限
            for(auto &U : availableU) lower = min(lower, distance(targetU, U) * distance(targetU, U));
            // std::cout << "2eps' " << 2*eps_prime << std::endl;
            // std::cout << "下限 " << lower << std::endl;
            // std::cout << "|X| " << availableU.size() << std::endl;
            if(availableU.empty()) return {};
            // if(availableU.size() > 120) std::cout << "|X| " << availableU.size() << std::endl;
            // for(auto U : availableU_ZOmega) std::cout << get_T_count(U) << " " << distance(targetU, to_quaternion(U)) << " " << to_quaternion(U) << std::endl;
            // FTYPE min_dist = 1.0;
            // for(int i = 0; i < availableU.size(); i++) for(int j = i+1; j < availableU.size(); j++) min_dist = min(min_dist, distance(availableU[i], availableU[j]));
            // std::cout << "min_dist " << min_dist << std::endl;

            if(target_eps < lower && !availableU.empty()) return {};   // 下限未満のときは、目的誤差は達成不可能なのでここで計算をやめる

            // std::cout << "eps_prime " << eps_prime << std::endl;
            // for(auto &U : availableU) std::cout << distance(targetU, U) << std::endl;

            // 1e-10を足しているのは数値誤差対策
            // std::cout << availableU.size() << std::endl;
            // if(eps_prime >= 1 || check_eps_net(availableU, targetU, eps_prime)){   // 2*eps_prime > 1のときは2ε'近傍がS^3全て
            if(2*eps_prime >= 0.99999 || availableU.size() >= 50){
                // std::cout << "bbb" << std::endl;
                auto [opt_dist, prob] = get_optimal_prob(availableU, targetU);
                // std::cout << "ccc" << std::endl;
                if(opt_dist > target_eps) return {};   // 最適な確率混合で目的誤差を達成できてないときは空を返す

                std::vector<std::pair<FTYPE, U2_ZOmega>> ret;
                for(int i = 0; i < prob.size(); i++){
                    if(prob[i] > 1e-15) ret.push_back({prob[i], availableU_ZOmega[i]});
                }
                return ret;
            }

            eps_prime += delta_eps_prime;
            eps_prime = min(eps_prime, (FTYPE)0.5);
            // std::cout << eps_prime << std::endl;
        }

    }




    std::vector<std::pair<FTYPE, std::string>> Prob_Unitary_Synthesis(quaternion targetU, FTYPE eps)
    {
        int MAX_T_count = 100;
        for(int T_count = 0; T_count < MAX_T_count; T_count++){
            // std::cout << std::endl;
            // std::cout << "Tカウント " << T_count << std::endl;
            // std::cout << std::endl;
            // auto optimal_mixed_unitary = optimal_prob_unitary(T_count, targetU, eps);
            auto optimal_mixed_unitary = optimal_prob_unitary(T_count, targetU, eps);
            if(!optimal_mixed_unitary.empty()){
                std::vector<std::pair<FTYPE, std::string>> ret;
                for(auto& [p, U] : optimal_mixed_unitary) ret.push_back({p, ExactSynthesis(U)});
                return ret;
            } 
        }

        return {};   // 警告出すとかが良いのか？
    }



    int get_T_count(std::vector<std::pair<FTYPE, std::string>> mixed_unitary)
    {
        int T_count = -1;
        for(auto [p, seq] : mixed_unitary) T_count = std::max(T_count, (int)std::count(seq.begin(), seq.end(), 'T'));
        return T_count;
    }
}

