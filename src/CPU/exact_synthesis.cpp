#include <bits/stdc++.h>
#include <Eigen/Core>
#include "Rings.cpp"
#include "quaternion.hpp"

// Bloch sphere representation
template <typename T>
std::pair<Eigen::Matrix<ZRoot2<T>, 3, 3>, int> to_SO3(Eigen::Matrix<ZOmega<T>, 2, 2> U, int k){
    Eigen::Matrix<ZOmega<T>, 2, 2> U_dag = adjoint(U);
    Eigen::Matrix<ZRoot2<T>, 3, 3> U_SO3;
    int k_SO3 = 2 * k + 1;

    ZOmega<T> zero = {0,0,0,0};
    ZOmega<T> one = {0,0,0,1};
    ZOmega<T> i_u = {0,1,0,0};
    ZOmega<T> sqrt2 = {-1,0,1,0};


    Eigen::Matrix<ZOmega<T>, 2, 2> X, Y, Z;
    X << zero, one, 
         one , zero;
    Y << zero, -i_u,
         i_u, zero;
    Z << one , zero,
         zero, -one;

    Eigen::Matrix<ZOmega<T>, 2, 2> UXU_dag, UYU_dag, UZU_dag;
    UXU_dag = sqrt2 * U * X * U_dag;
    UYU_dag = sqrt2 * U * Y * U_dag;
    UZU_dag = sqrt2 * U * Z * U_dag;

    U_SO3(0,0) = UXU_dag(0,1).real();
    U_SO3(0,1) = UYU_dag(0,1).real();
    U_SO3(0,2) = UZU_dag(0,1).real();

    U_SO3(1,0) = -UXU_dag(0,1).imag();
    U_SO3(1,1) = -UYU_dag(0,1).imag();
    U_SO3(1,2) = -UZU_dag(0,1).imag();

    U_SO3(2,0) = UXU_dag(0,0).real();
    U_SO3(2,1) = UYU_dag(0,0).real();
    U_SO3(2,2) = UZU_dag(0,0).real();
    
    reduction(U_SO3, k_SO3);

    return {U_SO3, k_SO3};
}

template <typename T, int N, int M>
Eigen::Matrix<int, N, M> parity(Eigen::Matrix<ZRoot2<T>, N, M> A){
    Eigen::Matrix<int, N, M> Ap;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++) Ap(i,j) = abs(A(i,j).a % 2);
    }
    return Ap;
}


std::vector<std::pair<std::string, Eigen::Matrix<int, 3, 3>>> Clifford;

void Generate_Clifford(){
    if(Clifford.size() == 24) return;
    Clifford.resize(24);

    Clifford[0].first = "";
    Clifford[1].second << 1, 0, 0,
                     0, 1, 0,
                     0, 0, 1;
    Clifford[1].first = "S";
    Clifford[1].second << 0, -1, 0, 
                     1, 0, 0, 
                     0, 0, 1;
    Clifford[2].first = "H";
    Clifford[2].second << 0, 0, 1, 
                     0, -1, 0, 
                     1, 0, 0;
    Clifford[3].first = "SS";
    Clifford[3].second << -1, 0, 0, 
                     0, -1, 0, 
                     0, 0, 1;
    Clifford[4].first = "HS";
    Clifford[4].second << 0, 0, 1, 
                     -1, 0, 0, 
                     0, -1, 0;
    Clifford[5].first = "SH";
    Clifford[5].second << 0, 1, 0, 
                     0, 0, 1, 
                     1, 0, 0;
    Clifford[6].first = "SSS";
    Clifford[6].second << 0, 1, 0, 
                     -1, 0, 0, 
                     0, 0, 1;
    Clifford[7].first = "HSS";
    Clifford[7].second << 0, 0, 1, 
                     0, 1, 0, 
                     -1, 0, 0;
    Clifford[8].first = "SHS";
    Clifford[8].second << 1, 0, 0, 
                     0, 0, 1, 
                     0, -1, 0;
    Clifford[9].first = "SSH";
    Clifford[9].second << 0, 0, -1, 
                     0, 1, 0, 
                     1, 0, 0;
    Clifford[10].first = "HSH";
    Clifford[10].second << 1, 0, 0, 
                     0, 0, -1, 
                     0, 1, 0;
    Clifford[11].first = "HSSS";
    Clifford[11].second << 0, 0, 1, 
                      1, 0, 0, 
                      0, 1, 0;
    Clifford[12].first = "SHSS";
    Clifford[12].second << 0, -1, 0, 
                      0, 0, 1, 
                      -1, 0, 0;
    Clifford[13].first = "SSHS";
    Clifford[13].second << 0, 0, -1, 
                      1, 0, 0, 
                      0, -1, 0;
    Clifford[14].first = "HSHS";
    Clifford[14].second << 0, -1, 0, 
                      0, 0, -1, 
                      1, 0, 0;
    Clifford[15].first = "HSSH";
    Clifford[15].second << 1, 0, 0, 
                      0, -1, 0, 
                      0, 0, -1;
    Clifford[16].first = "SHSSS";
    Clifford[16].second << -1, 0, 0, 
                      0, 0, 1, 
                      0, 1, 0;
    Clifford[17].first = "SSHSS";
    Clifford[17].second << 0, 0, -1, 
                      0, -1, 0, 
                      -1, 0, 0;
    Clifford[18].first = "HSHSS";
    Clifford[18].second << -1, 0, 0, 
                      0, 0, -1, 
                      0, -1, 0;
    Clifford[19].first = "HSSHS";
    Clifford[19].second << 0, -1, 0, 
                      -1, 0, 0, 
                      0, 0, -1;
    Clifford[20].first = "SHSSH";
    Clifford[20].second << 0, 1, 0, 
                      1, 0, 0, 
                      0, 0, -1;
    Clifford[21].first = "SSHSSS";
    Clifford[21].second << 0, 0, -1, 
                      -1, 0, 0, 
                      0, 1, 0;
    Clifford[22].first = "HSHSSS";
    Clifford[22].second << 0, 1, 0, 
                      0, 0, -1, 
                      -1, 0, 0;
    Clifford[23].first = "HSSHSS";
    Clifford[23].second << -1, 0, 0, 
                      0, 1, 0, 
                      0, 0, -1;
}


template <typename T>
std::string Exact_synthesis(ZOmega<T> u, ZOmega<T> t, int k){
    Eigen::Matrix<ZOmega<T>, 2, 2> U;
    U << u, -t.conj(), 
         t,  u.conj();
    
    Generate_Clifford();
    
    auto [U_SO3, k_SO3] = to_SO3(U, k);

    Eigen::Matrix<ZRoot2<T>, 3, 3> H_SO3, S_SO3, T_SO3;
    H_SO3(0,0) = H_SO3(0,1) = H_SO3(1,0) = H_SO3(1,2) = H_SO3(2,1) = H_SO3(2,2) = 0;
    H_SO3(0,2) = H_SO3(2,0) = 1;
    H_SO3(1,1) = -1;
    S_SO3(0,0) = S_SO3(0,2) = S_SO3(1,1) = S_SO3(1,2) = S_SO3(2,0) = S_SO3(2,1) = 0;
    S_SO3(1,0) = S_SO3(2,2) = 1;
    S_SO3(0,1) = -1;
    T_SO3(0,2) = T_SO3(1,2) = T_SO3(2,0) = T_SO3(2,1) = 0;
    T_SO3(0,0) = T_SO3(1,0) = T_SO3(1,1) = 1;
    T_SO3(0,1) = -1;
    T_SO3(2,2) = {0,1};

    Eigen::Matrix<ZRoot2<T>, 3, 3> Hinv_SO3, Sinv_SO3, Tinv_SO3;
    Hinv_SO3 = H_SO3;
    Sinv_SO3 = S_SO3.transpose();
    Tinv_SO3 = T_SO3.transpose();
    

    std::string ret = "";
    while(k_SO3 > 0){
        Eigen::Matrix<int, 3, 3> P = parity(U_SO3);
        if(P(2,0) == 0 && P(2,1) == 0 && P(2,2) == 0){
            ret += "T";
            U_SO3 = Tinv_SO3 * U_SO3;
        }else if(P(0,0) == 0 && P(0,1) == 0 && P(0,2) == 0){
            ret += "HT";
            U_SO3 = Tinv_SO3 * Hinv_SO3 * U_SO3;
        }else{
            ret += "SHT";
            U_SO3 = Tinv_SO3 * Hinv_SO3 * Sinv_SO3 * U_SO3;
        }
        k_SO3++;   // U_SO3にT_invを作用させているからインクリメント
        
        reduction(U_SO3, k_SO3);
    }

    for(auto [seq, C] : Clifford){
        bool f = true;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++) if(U_SO3(i,j) != C(i,j)) f = false;
        }
        
        if(f) ret += seq;
    }

    return ret;
}