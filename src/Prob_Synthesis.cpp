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
#include "enum_u_t.hpp"
#include "ExactSynthesis.hpp"
#include "eps_net_verification.hpp"
#include "sdpa_dd.hpp"


namespace SU2_Compiler
{
    template <UINT dim>
    Eigen::Matrix<CTYPE, dim*dim, dim*dim> Choi_Jamiolkowski(Eigen::Matrix<CTYPE, dim, dim> U)
    {
        Eigen::Matrix<CTYPE, dim, dim> U_dag;
        U_dag = U.adjoint();
        Eigen::Matrix<CTYPE, dim*dim, dim*dim> CJ_U;
        
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                CJ_U(Eigen::seq(i*dim, (i+1)*dim-1), Eigen::seq(j*dim, (j+1)*dim-1)) = U.col(i) * U_dag.row(j);
            }
        }

        return CJ_U;
    }

    
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
        C[0] = C[1] = C[2] = C[3] = 0.5;

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


    std::pair<std::vector< std::pair<FTYPE, U2_ZOmega> >, std::vector<U2_ZOmega>>
    optimal_prob_unitary(int T_count, quaternion targetU, FTYPE eps, std::vector<U2_ZOmega> pre_availableU_ZOmega)
    {
        int k = (T_count + 1) / 2 + 1;

        FTYPE eps_prime = pow(2.0, -(FTYPE)T_count / 3.0) / 2.0;   // T_count = 3log_2(1/2eps')をeps'について解いた
        FTYPE delta_eps_prime = eps_prime / 5;

        std::vector<U2_ZOmega> new_availableU;
        std::array<std::vector<ZRoot2>, 12> pre_resluts = {};

        while(true){
            std::tie(new_availableU, pre_resluts) = enum_u_t(targetU, 2.0*eps_prime, k, T_count % 2, pre_resluts);
            
            std::vector<quaternion> availableU;
            std::vector<U2_ZOmega> availableU_ZOmega;
            for(auto U_ZOmega : new_availableU){
                if(get_T_count(U_ZOmega) > T_count) continue; 
                
                quaternion U = to_quaternion(U_ZOmega);
                availableU.push_back(U);
                availableU_ZOmega.push_back(U_ZOmega);
            }

            // T_count-1の結果なので、偶奇が逆になる
            for(auto U_ZOmega : pre_availableU_ZOmega){
                quaternion U = to_quaternion(U_ZOmega);

                if(distance(targetU, U) > 2*eps_prime) continue; 

                availableU.push_back(U);
                availableU_ZOmega.push_back(U_ZOmega);
            }
            
            std::cout << "eps' " << eps_prime << std::endl;
            std::cout << "利用可能なユニタリの数" << availableU.size() << std::endl;
            // std::sort(availableU.begin(), availableU.end(), [](quaternion const& lhs, quaternion const& rhs) {return lhs.a < rhs.a;});
            // for(auto x : availableU) std::cout << std::setprecision(30) << x << std::endl;
            // exit(0);

            FTYPE min_dist = 1.0;
            for(auto u : availableU) min_dist = std::min(min_dist, distance(targetU, u));
            // if(!availableU.empty() && min_dist * min_dist > eps) return {{}, {}};

            if(eps_prime > 1 || check_eps_net(availableU, targetU, eps_prime)){
                auto [opt_dist, prob] = get_optimal_prob(availableU, targetU);
                for(auto p : prob) std::cout << p << " ";
                std::cout << std::endl;
                std::cout << "下限 " << min_dist*min_dist << std::endl;
                // std::cout << std::setprecision(20) << diamond_distance(targetU, availableU, prob) << std::endl;

                std::vector<std::pair<FTYPE, U2_ZOmega>> ret;
                for(int i = 0; i < prob.size(); i++){
                    if(prob[i] > 1e-20) ret.push_back({prob[i], availableU_ZOmega[i]});
                }
                return {ret, new_availableU};
            }
            
            eps_prime += delta_eps_prime;
        }
    }



    std::vector<std::pair<FTYPE, std::string>> Prob_Unitary_Synthesis(quaternion targetU, FTYPE eps)
    {
        FTYPE eps_sqrt = sqrt(eps);

        int T_count = 0;
        std::vector< std::pair<FTYPE, U2_ZOmega> > optimal_CliffordT;
        std::vector<U2_ZOmega> pre_availableU_ZOmega = {};
        while(true){
            std::tie(optimal_CliffordT, pre_availableU_ZOmega) = optimal_prob_unitary(T_count, targetU, eps, pre_availableU_ZOmega);
            
            if(distance(targetU, optimal_CliffordT) < eps)
                
            T_count++;
        }
    }

}

