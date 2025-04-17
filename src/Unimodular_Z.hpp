#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "type.hpp"
#include "rings.hpp"
#include "grid_solver.hpp"

namespace SU2_Compiler
{

template <int N>
using UnimodularMatrix = Eigen::Matrix<ITYPE, N, N>;


template <int N>
std::ostream& operator<<(std::ostream& os, const UnimodularMatrix<N>& U){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) os << U(i,j) << " ";
        if(i != N-1) os << std::endl;
    }
    return os;
}


template <typename T>
std::vector<std::vector<T>> getSubMatrix(std::vector<std::vector<T>>& A, int p, int q, int N){
    std::vector<std::vector<T>> A_sub(N-1, std::vector<T>(N-1));
    int i = 0, j = 0;
    for(int row = 0; row < N; row++){
        for(int col = 0; col < N; col++){
            if(row == p || col == q) continue;

            A_sub[i][j] = A[row][col];
            j++;
            if(j == N-1){
                j = 0;
                i++;
            }
        }
    }
    return A_sub;
}

template <typename T>
T determinant(std::vector<std::vector<T>> A, int N){
    if(N == 1) return A[0][0];
    
    T ret = 0;
    for(int i = 0; i < N; i++){
        std::vector<std::vector<T>> A_sub = getSubMatrix(A, 0, i, N);
        if(i % 2) ret -= A[0][i] * determinant(A_sub, N-1);
        else ret += A[0][i] * determinant(A_sub, N-1);
    }
    return ret;
} 

template <typename T, int N>
T determinant(Eigen::Matrix<T, N, N> A){
    std::vector<std::vector<T>> A_vec(N, std::vector<T>(N));
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) A_vec[i][j] = A(i,j);
    }

    return determinant(A_vec, N);
}


// ユニモジュラ行列の逆行列計算
template <int N>
UnimodularMatrix<N> inverse_UnimodularMatrix(UnimodularMatrix<N>& U){
    std::vector<std::vector<ITYPE>> U_vec(N, std::vector<ITYPE>(N));
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) U_vec[i][j] = U(i,j);
    }
    
    UnimodularMatrix<N> U_inv;
    ITYPE detU = determinant(U_vec, N);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            auto subU_ij = getSubMatrix<ITYPE>(U_vec, i, j, N);
            if((i + j) % 2) U_inv(j,i) = -detU * determinant(subU_ij, N-1);
            else U_inv(j,i) = detU * determinant(subU_ij, N-1);
        }
    }
    return U_inv;
}


template <int M, int N>
Eigen::Matrix<FTYPE, M, N> GSO(const Eigen::Matrix<FTYPE, M, N>& B){
    Eigen::Matrix<FTYPE, M, N> B_aste;
    auto mu = [&](int i, int j) { return B.col(i).dot(B_aste.col(j)) / B_aste.col(j).dot(B_aste.col(j)); };

    for(int i = 0; i < N; i++){
        B_aste.col(i) = B.col(i);
        for(int j = 0; j < i; j++){
            B_aste.col(i) -= mu(i, j) * B_aste.col(j);
        }
    }

    return B_aste;
} 


#include <chrono>

template <int M, int N>
UnimodularMatrix<N> LLL(Eigen::Matrix<FTYPE, M, N> B, FTYPE delta){
    using namespace std;
    chrono::system_clock::time_point start, end;
    double time;

    Eigen::Matrix<FTYPE, M, N> B_aste;
    Eigen::Matrix<FTYPE, N, N> mu;
    Eigen::Vector<FTYPE, N> abs_B_aste;
    B_aste = GSO<M, N>(B);
    for(int i = 0; i < N; i++) abs_B_aste(i) = B_aste.col(i).dot(B_aste.col(i));
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) mu(i,j) = B.col(i).dot(B_aste.col(j)) / B_aste.col(j).dot(B_aste.col(j));
    }

    UnimodularMatrix<N> U;
    U = UnimodularMatrix<N>::Identity();

    int cnt = 0;
    int k = 1;
    while(k < N){
        // std::cout << "k = " << k << std::endl;
        cnt++;
        if(cnt >= 100000) std::cout << "LLLが無限ループ" << std::endl;

        start = chrono::system_clock::now();
        // Size reduce operation
        for(int j = k-1; j >= 0; j--){
            if(abs(mu(k,j)) > 0.5){
                ITYPE q = (ITYPE)round(mu(k,j));
                B.col(k) -= q * B.col(j);
                U.col(k) -= q * U.col(j);
                for(int l = 0; l <= j; l++) mu(k,l) -= q * mu(j,l);
                // B_aste = GSO(B);
                // for(int i = 0; i < N; i++) abs_B_aste(i) = B_aste.col(i).dot(B_aste.col(i));
                // for(int i = 0; i < N; i++){
                //     for(int j = 0; j < N; j++) mu(i,j) = B.col(i).dot(B_aste.col(j)) / B_aste.col(j).dot(B_aste.col(j));
                // }

            }
        }


        if(abs_B_aste(k) >= (delta - mu(k,k-1)*mu(k,k-1)) * abs_B_aste(k-1)){
            k++;
        }else{
            B.col(k-1).swap(B.col(k));
            U.col(k-1).swap(U.col(k));

            FTYPE mu_prime = mu(k,k-1);
            FTYPE new_abs_B = abs_B_aste(k) + mu_prime*mu_prime * abs_B_aste(k-1);
            mu(k,k-1) = mu_prime * abs_B_aste(k-1) / new_abs_B;
            abs_B_aste(k) = abs_B_aste(k) * abs_B_aste(k-1) / new_abs_B;
            abs_B_aste(k-1) = new_abs_B;

            for(int j = 0; j <= k-2; j++){
                std::swap(mu(k-1,j), mu(k,j));
            }

            for(int j = k+1; j < N; j++){
                FTYPE t = mu(j,k);
                mu(j,k) = mu(j,k-1) - mu_prime * t;
                mu(j,k-1) = t + mu(k,k-1) * mu(j,k);
            }

            // B_aste = GSO(B);
            // for(int i = 0; i < N; i++) abs_B_aste(i) = B_aste.col(i).dot(B_aste.col(i));
            // for(int i = 0; i < N; i++){
            //     for(int j = 0; j < N; j++) mu(i,j) = B.col(i).dot(B_aste.col(j)) / B_aste.col(j).dot(B_aste.col(j));
            // }

            k = std::max(k-1, 1);

            // FTYPE D = 1.0;
            // for(int i = 0; i < N-1; i++){
            //     FTYPE d_i = 1.0;
            //     for(int j = 0; j <= i; j++) d_i *= abs_B_aste(j);
            //     D *= d_i;
            // }
            // std::cout << D << std::endl;
        }
    }

    // std::cout << "X_reduced = \n" << X << std::endl;
    // FTYPE r = 1.0;
    // for(int i = 0; i < N; i++){
    //     r *= B.col(i).dot(B.col(i)) / abs_B_aste(i);
    // }
    // std::cout << r << std::endl;

    // std::cout << "ループ数 " << cnt << std::endl;

    return U;
}



}