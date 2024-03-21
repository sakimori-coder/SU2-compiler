#include <bits/stdc++.h>
#include <Eigen/Core>
#include <boost/math/constants/constants.hpp>
#include "quaternion.hpp"
#include "sdpa_dd.hpp"
#include "enumerate_u_t.cpp"
#include "eps_net_verification.cpp"


using Unitary_ZOmega = std::tuple<ZOmega<long long>, ZOmega<long long>, int>;   // u, t, kのタプル


// UnitaryのChoi行列
template <typename T, int dim>
Eigen::Matrix<T, dim*dim, dim*dim> Choi_Jamiolkowski(Eigen::Matrix<T, dim, dim> U){
    Eigen::Matrix<T, dim, dim> U_dag;
    U_dag = U.adjoint();
    Eigen::Matrix<T, dim*dim, dim*dim> CJ_U;
    
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            CJ_U(Eigen::seq(i*dim, (i+1)*dim-1), Eigen::seq(j*dim, (j+1)*dim-1)) = U.col(i) * U_dag.row(j);
        }
    }

    return CJ_U;
}


// UnitaryのChoi行列のMagic基底表現
template <typename T>
Eigen::Matrix<T, 4, 4> CJ_MB(Eigen::Matrix<std::complex<T>, 2, 2> U){
    Eigen::Matrix<std::complex<T>, 4, 4> CJ_U = Choi_Jamiolkowski<std::complex<T>, 2>(U);

    std::complex<T> i(0., 1.);
    T sqrt2_inv = 1.0 / sqrt((T)2.0);
    std::complex<T> isqrt2_inv(0, sqrt2_inv);
    Eigen::Matrix<std::complex<T>, 4, 4> M = Eigen::Matrix<std::complex<T>, 4, 4>::Zero(4,4);
    M(0,0) = sqrt2_inv;
    M(0,3) = sqrt2_inv;
    M(1,0) =-isqrt2_inv;
    M(1,3) = isqrt2_inv;
    M(2,1) =-isqrt2_inv;
    M(2,2) =-isqrt2_inv;
    M(3,1) = sqrt2_inv;
    M(3,2) =-sqrt2_inv;

    Eigen::Matrix<std::complex<T>, 4, 4> CJ_U_MB;
    CJ_U_MB = M * CJ_U * M.adjoint();

    return CJ_U_MB.real();
}


template <typename T, int N, int M>
std::vector<std::vector<T>> to_2d_vec(Eigen::Matrix<T, N, M> A){
    std::vector<std::vector<T>> ret(N, std::vector<T>(M));
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++) ret[i][j] = A(i,j);
    }
    return ret;
}


// ユニタリと確率混合ユニタリのダイアモンド距離を計算
template <typename T>
T diamond_distance(quaternion<T> A, std::vector<quaternion<T>> B, std::vector<T> prob){
    const int N = B.size();
    
    Eigen::Matrix<T, 4, 4> CJ_A_MB = CJ_MB<T>(A.get_Matrix());
    std::vector<Eigen::Matrix<T, 4, 4>> CJ_B_MB(N);
    for(int i = 0; i < N; i++) CJ_B_MB[i] = CJ_MB<T>(B[i].get_Matrix());

    std::vector<Eigen::Matrix<T, 4, 4>> Symmetric_basis(10, Eigen::Matrix<T, 4, 4>::Zero(4,4));
    for(int i = 0; i < 4; i++) Symmetric_basis[i](i,i) = 1;
    Symmetric_basis[4](0,1) = Symmetric_basis[4](1,0) = 1;
    Symmetric_basis[5](0,2) = Symmetric_basis[5](2,0) = 1;
    Symmetric_basis[6](0,3) = Symmetric_basis[6](3,0) = 1;
    Symmetric_basis[7](1,2) = Symmetric_basis[7](2,1) = 1;
    Symmetric_basis[8](1,3) = Symmetric_basis[8](3,1) = 1;
    Symmetric_basis[9](2,3) = Symmetric_basis[9](3,2) = 1;

    int num_variable = 10;

    std::vector<T> C(num_variable, 0);
    C[0] = C[1] = C[2] = C[3] = 0.5;

    std::vector<std::vector<T>> ZERO_mat(4, std::vector<T>(4, 0));
    vector_4d<T> F(num_variable + 1, std::vector<std::vector<std::vector<T>>>(2));
    
    F[0][0] = ZERO_mat;
    Eigen::Matrix<T, 4, 4> D = -CJ_A_MB;
    for(int i = 0; i < N; i++) D += prob[i]*CJ_B_MB[i];
    F[0][1] = to_2d_vec<T, 4, 4>(D);

    for(int i = 1; i <= 10; i++){
        F[i][0] = to_2d_vec<T, 4, 4>(Symmetric_basis[i-1]);
        F[i][1] = to_2d_vec<T, 4, 4>(Symmetric_basis[i-1]);
    } 

    std::vector<T> xVec = sdpa_solver(C, F);

    return (xVec[0] + xVec[1] + xVec[2] + xVec[3]) / 2;
}





template <typename T>
std::vector<T> get_optimal_prob(std::vector<quaternion<T>> availableU, quaternion<T> targetU){
    const int N = availableU.size();
    
    std::vector<Eigen::Matrix<T, 4, 4>> availableCJUMB(N);
    for(int i = 0; i < N; i++) availableCJUMB[i] = CJ_MB<T>(availableU[i].get_Matrix());
    Eigen::Matrix<T, 4, 4> targetCJUMB = CJ_MB<T>(targetU.get_Matrix());

    
    std::vector<Eigen::Matrix<T, 4, 4>> Symmetric_basis(10, Eigen::Matrix<T, 4, 4>::Zero(4,4));
    for(int i = 0; i < 4; i++) Symmetric_basis[i](i,i) = 1;
    Symmetric_basis[4](0,1) = Symmetric_basis[4](1,0) = 1;
    Symmetric_basis[5](0,2) = Symmetric_basis[5](2,0) = 1;
    Symmetric_basis[6](0,3) = Symmetric_basis[6](3,0) = 1;
    Symmetric_basis[7](1,2) = Symmetric_basis[7](2,1) = 1;
    Symmetric_basis[8](1,3) = Symmetric_basis[8](3,1) = 1;
    Symmetric_basis[9](2,3) = Symmetric_basis[9](3,2) = 1;
    
    int num_variable = 10 + N;

    // 変数の並び順は (Sの変数10個, 混合確率)
    std::vector<T> C(num_variable, 0);
    C[0] = C[1] = C[2] = C[3] = 0.5;

    std::vector<std::vector<T>> ZERO_mat(4, std::vector<T>(4, 0));
    std::vector<std::vector<T>> ZERO_vec(1, std::vector<T>(N, 0)); 
    vector_4d<T> F(num_variable + 1, std::vector<std::vector<std::vector<T>>>(4));
    
    // F0
    F[0][0] = ZERO_mat;
    F[0][1] = to_2d_vec<T, 4, 4>(targetCJUMB);
    F[0][2] = ZERO_vec;
    F[0][3] = {{-1}};

    // F1~F10
    for(int i = 1; i <= 10; i++){
        F[i][0] = to_2d_vec<T, 4, 4>(Symmetric_basis[i-1]);
        F[i][1] = to_2d_vec<T, 4, 4>(Symmetric_basis[i-1]);
        F[i][2] = ZERO_vec;
        F[i][3] = {{0}};
    }

    // F11~FN
    for(int i = 0; i < N; i++){
        F[i+11][0] = ZERO_mat;
        F[i+11][1] = to_2d_vec<T, 4, 4>(availableCJUMB[i]);
        F[i+11][2] = ZERO_vec; 
        F[i+11][2][0][i] = 1;
        F[i+11][3] = {{-1}};
    }
    
    std::vector<T> xVec = sdpa_solver(C, F);    
    std::vector<T> prob(N);
    for(int i = 0; i < N; i++) prob[i] = xVec[i+10];

    T opt_val = (xVec[0] + xVec[1] + xVec[2] + xVec[3]) / 2;

    return prob;
}



template <typename T>
std::pair<std::vector<std::pair<T, Unitary_ZOmega>>, std::vector<Unitary_ZOmega>>
optimal_prob_unitary(int T_count, quaternion<T> targetU, std::vector<Unitary_ZOmega> pre_availableU = {})
{
    int k = (T_count + 1) / 2 + 1;

    T eps = pow(2.0, -(T)T_count / 3.0) / 2.0;   // T_count = 3log_2(1/2eps)をepsについて解いた
    T delta_eps = eps / 5;

    quaternion<T> omega_sqrt(cos(boost::math::constants::pi<T>() / 8), sin(boost::math::constants::pi<T>() / 8), 0, 0);
    std::vector<Unitary_ZOmega> new_availableU;
    std::array<std::vector<ZRoot2<long long>>, 12> pre_resluts = {};

    while(true){
        std::tie(new_availableU, pre_resluts) = enumerate_u_t_wrapper<long long, T>(targetU, 2.0*eps, k, pre_resluts, T_count % 2);
        
        std::vector<quaternion<T>> availableU;
        std::vector<Unitary_ZOmega> availableU_ZOmega;
        for(auto [u, t, k] : new_availableU){
            // std::string str = Exact_synthesis(u, t, k, T_count % 2);
            if(get_T_count(u, t, k, T_count % 2) > T_count) continue; 

            // cout << std::count(str.begin(), str.end(), 'T') << endl;
            // cout << u << endl;
            // cout << t << endl;
            // cout << distance(targetU, to_quaternion<T>(str)) << endl;
            
            quaternion<T> U = to_quaterion<long long, T>(u, t);
            U = U / U.norm();
            if(T_count % 2) U = omega_sqrt.conj() * U;
            availableU.push_back(U);
            availableU_ZOmega.push_back({u, t, k});
        }

        // T_count-1の結果なので、偶奇が逆になる
        for(auto [u, t, pre_k] : pre_availableU){
            if(get_T_count(u, t, k, (T_count-1) % 2) > T_count) continue;   // このif文には入らないと思う

            quaternion<T> U = to_quaterion<long long, T>(u, t);
            U = U / U.norm();
            if((T_count-1) % 2) U = omega_sqrt.conj() * U;
            availableU.push_back(U);
            availableU_ZOmega.push_back({u, t, pre_k});
        }
        
        std::cout << "eps " << eps << std::endl;
        std::cout << "利用可能なユニタリの数" << availableU.size() << std::endl;
        std::sort(availableU.begin(), availableU.end(), [](quaternion<T> const& lhs, quaternion<T> const& rhs) {return lhs.get_a() < rhs.get_a();});
        // for(auto x : availableU) std::cout << std::setprecision(30) << x << std::endl;
        // exit(0);

        if(check_eps_net(availableU, targetU, eps)){
            T min_eps = 1.0;
            for(auto u : availableU) min_eps = std::min(min_eps, distance(targetU, u));
            std::vector<T> prob = get_optimal_prob<T>(availableU, targetU);
            for(auto p : prob) std::cout << p << " ";
            std::cout << std::endl;
            std::cout << "下限 " << min_eps*min_eps << std::endl;
            // std::cout << std::setprecision(20) << diamond_distance(targetU, availableU, prob) << std::endl;

            std::vector<std::pair<T, Unitary_ZOmega>> ret;
            for(int i = 0; i < prob.size(); i++){
                if(prob[i] > 1e-20) ret.push_back({prob[i], availableU_ZOmega[i]});
            }
            return {ret, new_availableU};
        }
        
        eps += delta_eps;
    }
}


template <typename T>
std::vector<T, std::string> Prob_unitary_synthesis(quaternion<T> targetU, T eps){
    T eps_sqrt = sqrt(eps);
    T flag_lower_bound = false;

    int T_count = 0;
    std::vector<std::pair<T, Unitary_ZOmega> ans;
    std::vector<Unitary_ZOmega> pre_availableU = {};
    while(true){
        // 下限がターゲットユニタリに届いていなかった場合
        if(!flag_lower_bound){
            quaternion<T> targetU_prime = targetU;
            int k = (T_count + 1) / 2 + 1;
            auto [availableU, garbage] = enumerate_u_t<long long, T>(targetU, eps_sqrt, k, pre_resluts, T_count % 2);
            if(!availableU.empty()){
                flag_lower_bound = true;
                std::tie(ans, pre_availableU) = optimal_prob_unitary(T_count, targetU);
            } 
        }

        if(flag_lower_bound){
            std::tie(ans, pre_availableU) = optimal_prob_unitary(T_count, targetU, pre_availableU);

            for(auto [u, t, k] : ans){
                if(k == (T_count + 1) / 2 + 1){
                    quaternion<T> U = to_quaterion<long long, T>(u, t);
                    U = U / U.norm();
                    if(T_count % 2) U = omega_sqrt.conj() * U;
                    
                }else{
                    quaternion<T> U = to_quaterion<long long, T>(u, t);
                    U = U / U.norm();
                    if((T_count-1) % 2) U = omega_sqrt.conj() * U;              
                }
            }
        }
    }
}