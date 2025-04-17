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
using UnimodularMatrix = Eigen::Matrix<ZRoot2, N, N>;

template <int N>
UnimodularMatrix<N> conj(UnimodularMatrix<N> U){
    UnimodularMatrix<N> U_dot;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) U_dot(i,j) = conj(U(i,j));
    }
    return U_dot;
}

template <int N>
Eigen::Matrix<FTYPE, N, N> ZRoot2_to_FTYPE(UnimodularMatrix<N> U){
    Eigen::Matrix<FTYPE, N, N> U_FTYPE;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) U_FTYPE(i,j) = ZRoot2_to_FTYPE(U(i,j));
    }
    return U_FTYPE;
}

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
    std::vector<std::vector<ZRoot2>> U_vec(N, std::vector<ZRoot2>(N));
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) U_vec[i][j] = U(i,j);
    }
    
    UnimodularMatrix<N> U_inv;
    ZRoot2 detU = determinant(U_vec, N);
    ZRoot2 inv_detU = -conj(detU);
    if(detU * inv_detU == ZRoot2(-1,0)) inv_detU = -inv_detU;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            auto subU_ij = getSubMatrix<ZRoot2>(U_vec, i, j, N);
            if((i + j) % 2) U_inv(j,i) = -inv_detU * determinant(subU_ij, N-1);
            else U_inv(j,i) = inv_detU * determinant(subU_ij, N-1);
        }
    }
    return U_inv;
}


template <int N>
Eigen::Matrix<FTYPE, N, N> GSO(const Eigen::Matrix<FTYPE, N, N>& B){
    Eigen::Matrix<FTYPE, N, N> B_aste;
    auto mu = [&](int i, int j) { return B.col(i).dot(B_aste.col(j)) / B_aste.col(j).dot(B_aste.col(j)); };

    for(int i = 0; i < N; i++){
        B_aste.col(i) = B.col(i);
        for(int j = 0; j < i; j++){
            B_aste.col(i) -= mu(i, j) * B_aste.col(j);
        }
    }

    return B_aste;
} 


FTYPE pow_int(FTYPE x, int n){
    if(n < 0){
        n = -n;
        x = 1.0/x;
    }
    FTYPE ret = 1.0;
    while(n > 0){
        if(n & 1) ret *= x;
        x *= x;
        n >>= 1;
    }
    return ret;
}


std::pair<ITYPE, ITYPE> CVP_two_dimension(
    Eigen::Vector<FTYPE, 2> b1, 
    Eigen::Vector<FTYPE, 2> b2, 
    Eigen::Vector<FTYPE, 2> t)
{
    Eigen::Matrix<ITYPE, 2, 2> U;
    U = Eigen::Matrix<ITYPE, 2, 2>::Identity();
    FTYPE B1 = b1.dot(b1);
    FTYPE B2 = b2.dot(b2);
    while(true){
        if(B2 < B1){
            std::swap(b1, b2);
            U.col(0).swap(U.col(1));
            B1 = B2;
        }

        ITYPE mu = (ITYPE)round(b1.dot(b2) / B1);
        b2 -= (FTYPE)mu * b1;
        B2 = b2.dot(b2);
        U.col(1) += mu * U.col(0);
        if(B2 >= B1) break;
    }

    Eigen::Vector<ITYPE, 2> ab;
    Eigen::Vector<FTYPE, 2> b2_aste = b2 - b2.dot(b1) / B1 * b1;
    ab(1) = (ITYPE)round(b2_aste.dot(t) / b2_aste.dot(b2_aste));
    ab(0) = (ITYPE)round(b1.dot(t - (FTYPE)ab(1) * b2) / B1);
    ITYPE det = U(0,0) * U(1,1) - U(0,1) * U(1,0);
    Eigen::Matrix<ITYPE, 2, 2> invU;
    invU << U(1,1),-U(0,1),
           -U(1,0), U(0,0);
    invU *= det;
    ab = invU * ab;
    return {ab(0), ab(1)};
}

#include <chrono>

template <int N>
UnimodularMatrix<N> LLL_ZRoot2(Eigen::Matrix<FTYPE, N, N> X, Eigen::Matrix<FTYPE, N, N> Y, FTYPE delta){
    using namespace std;
    chrono::system_clock::time_point start, end;
    double time;
    // std::cout << "X = \n" << X << std::endl;
    // std::cout << "Y = \n" << Y << std::endl;
    Eigen::Matrix<FTYPE, N, N> X_aste, Y_aste;
    X_aste = GSO(X);
    Y_aste = GSO(Y);
    auto mu = [&](int i, int j) { return X.col(i).dot(X_aste.col(j)) / X_aste.col(j).dot(X_aste.col(j)); };
    auto nu = [&](int i, int j) { return Y.col(i).dot(Y_aste.col(j)) / Y_aste.col(j).dot(Y_aste.col(j)); };
    auto abs_X_aste = [&](int i) { return X_aste.col(i).dot(X_aste.col(i)); };
    auto abs_Y_aste = [&](int i) { return Y_aste.col(i).dot(Y_aste.col(i)); };

    FTYPE C = (sqrt2 / 2) * (sqrt2 / 2) + ((1 + sqrt2) / 2) * ((1 + sqrt2) / 2);
    UnimodularMatrix<N> U;
    U = UnimodularMatrix<N>::Identity();

    
    // Normalize操作は最初の一回だけはFTYPE型で実行. 今後はdouble型で実行する
    // FTYPE型のlogが重い
    for(int i = 0; i < N; i++){
        FTYPE r = abs_X_aste(i) / abs_Y_aste(i);
        int m = -round((log(r) / log(LAMBDA)) / 4.0);
        ZRoot2 LAMBDA_pow_m = pow_lambda(m);
        ZRoot2 LAMBDA_dot_pow_m = conj(LAMBDA_pow_m);
        X.col(i) *= ZRoot2_to_FTYPE(LAMBDA_pow_m);
        Y.col(i) *= ZRoot2_to_FTYPE(LAMBDA_dot_pow_m);
        X_aste.col(i) *= ZRoot2_to_FTYPE(LAMBDA_pow_m);
        Y_aste.col(i) *= ZRoot2_to_FTYPE(LAMBDA_dot_pow_m);
        U.col(i) *= LAMBDA_pow_m;
    }


    int cnt = 0;
    int k = 1;
    while(k < N){
        // std::cout << "k = " << k << std::endl;
        cnt++;
        if(cnt >= 100000) std::cout << "LLLが無限ループ" << std::endl;

        // start = chrono::system_clock::now();
        // Normalize operaton
        const double log_lambda = std::log((double)LAMBDA);
        for(int i = k-1; i <= k; i++){
            double r = abs_X_aste(i) / abs_Y_aste(i);
            int m = -std::round((std::log(r) / log_lambda) / 4.0);
            if(m == 0) continue;
            ZRoot2 LAMBDA_pow_m = pow_lambda(m);
            ZRoot2 LAMBDA_dot_pow_m = conj(LAMBDA_pow_m);
            X.col(i) *= ZRoot2_to_FTYPE(LAMBDA_pow_m);
            Y.col(i) *= ZRoot2_to_FTYPE(LAMBDA_dot_pow_m);
            X_aste.col(i) *= ZRoot2_to_FTYPE(LAMBDA_pow_m);
            Y_aste.col(i) *= ZRoot2_to_FTYPE(LAMBDA_dot_pow_m);
            U.col(i) *= LAMBDA_pow_m;
        }
        // end = chrono::system_clock::now();
        // time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
        // cout << "time1 " << time << "[ms]" << endl;


        start = chrono::system_clock::now();
        // Size reduce operation
        for(int j = k-1; j >= 0; j--){
            FTYPE mu_kj = mu(k, j);
            FTYPE nu_kj = nu(k, j);


            FTYPE c = sqrt(sqrt( abs_X_aste(j) / abs_Y_aste(j) ));
            if(c*c*mu_kj*mu_kj + nu_kj*nu_kj/(c*c) < 1.5) continue;
            Eigen::Vector<FTYPE, 2> b1, b2, t;
            b1 << c, 1.0/c;
            b2 << c*sqrt2, -1.0/c * sqrt2;
            t << c * (FTYPE)mu_kj, 1.0/c * (FTYPE)nu_kj;
            // std::cout << t << std::endl;
            // std::cout << std::endl;
            auto [a,b] = CVP_two_dimension(b1, b2, t);
            ZRoot2 mu_round(a,b);
            ZRoot2 nu_round = conj(mu_round);
            X.col(k) -= X.col(j) * ZRoot2_to_FTYPE(mu_round);
            Y.col(k) -= Y.col(j) * ZRoot2_to_FTYPE(nu_round);
            U.col(k) -= mu_round * U.col(j);
        }
        // end = chrono::system_clock::now();
        // time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
        // cout << "time2 " << time << "[ms]" << endl;


        // start = chrono::system_clock::now();
        FTYPE x_k_1 = abs_X_aste(k-1);
        FTYPE y_k_1 = abs_Y_aste(k-1);
        FTYPE x_k = abs_X_aste(k);
        FTYPE y_k = abs_Y_aste(k);
        FTYPE mu_kk_1 = mu(k,k-1);
        FTYPE nu_kk_1 = nu(k,k-1);
        if(delta*delta * x_k_1 * y_k_1 <= (x_k + mu_kk_1*mu_kk_1 * x_k_1) * (y_k + nu_kk_1*nu_kk_1 * y_k_1)){
            k++;
        }else{
            X.col(k-1).swap(X.col(k));
            Y.col(k-1).swap(Y.col(k));
            U.col(k-1).swap(U.col(k));
            Eigen::Vector<FTYPE, N> new_x_aste_k_1, new_x_aste_k, new_y_aste_k_1, new_y_aste_k;
            new_x_aste_k_1 = X_aste.col(k) + mu_kk_1 * X_aste.col(k-1);
            new_x_aste_k = x_k * X_aste.col(k-1) - mu_kk_1 * x_k_1 * X_aste.col(k);
            new_x_aste_k /= (x_k + mu_kk_1*mu_kk_1 * x_k_1);
            X_aste.col(k-1) = new_x_aste_k_1;
            X_aste.col(k) = new_x_aste_k;

            new_y_aste_k_1 = Y_aste.col(k) + nu_kk_1 * Y_aste.col(k-1);
            new_y_aste_k = y_k * Y_aste.col(k-1) - nu_kk_1 * y_k_1 * Y_aste.col(k);
            new_y_aste_k /= (y_k + nu_kk_1*nu_kk_1 * y_k_1);
            Y_aste.col(k-1) = new_y_aste_k_1;
            Y_aste.col(k) = new_y_aste_k;

            k = std::max(k-1, 1);
        }
        // end = chrono::system_clock::now();
        // time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
        // cout << "time3 " << time << "[ms]" << endl;
    }

    // std::cout << "X_reduced = \n" << X << std::endl;
    // std::cout << "Y_reduced = \n" << Y << std::endl;

    return U;
}



}