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
using GridOp = Eigen::Matrix<ZRoot2, N, N>;

template <int N>
GridOp<N> conj(GridOp<N> G){
    GridOp<N> G_dot;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) G_dot(i,j) = conj(G(i,j));
    }
    return G_dot;
}

template <int N>
Eigen::Matrix<FTYPE, N, N> ZRoot2_to_FTYPE(GridOp<N> G){
    Eigen::Matrix<FTYPE, N, N> G_FTYPE;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) G_FTYPE(i,j) = ZRoot2_to_FTYPE(G(i,j));
    }
    return G_FTYPE;
}

template <int N>
std::ostream& operator<<(std::ostream& os, const GridOp<N>& G){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) os << G(i,j) << " ";
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


// |det(G)|=1であるGridOpの逆行列計算
template <int N>
GridOp<N> inverse_Special_GridOp(GridOp<N>& G){
    std::vector<std::vector<ZRoot2>> G_vec(N, std::vector<ZRoot2>(N));
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) G_vec[i][j] = G(i,j);
    }
    
    GridOp<N> G_inv;
    ZRoot2 detG = determinant(G_vec, N);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            auto G_ij = getSubMatrix<ZRoot2>(G_vec, i, j, N);
            if((i + j) % 2) G_inv(j,i) = -detG * determinant(G_ij, N-1);
            else G_inv(j,i) = detG * determinant(G_ij, N-1);
        }
    }
    return G_inv;
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


template <int N>
GridOp<N> LLL_ZRoot2(Eigen::Matrix<FTYPE, N, N> X, Eigen::Matrix<FTYPE, N, N> Y, FTYPE delta){
    std::cout << "X = \n" << X << std::endl;
    std::cout << "Y = \n" << Y << std::endl;
    Eigen::Matrix<FTYPE, N, N> X_aste, Y_aste;
    X_aste = GSO(X);
    Y_aste = GSO(Y);
    auto mu = [&](int i, int j) { return X.col(i).dot(X_aste.col(j)) / X_aste.col(j).dot(X_aste.col(j)); };
    auto nu = [&](int i, int j) { return Y.col(i).dot(Y_aste.col(j)) / Y_aste.col(j).dot(Y_aste.col(j)); };
    auto abs_X_aste = [&](int i) { return X_aste.col(i).dot(X_aste.col(i)); };
    auto abs_Y_aste = [&](int i) { return Y_aste.col(i).dot(Y_aste.col(i)); };

    FTYPE C = (sqrt2 / 2) * (sqrt2 / 2) + ((1 + sqrt2) / 2) * ((1 + sqrt2) / 2);
    GridOp<N> G;
    G = GridOp<N>::Identity();

    int k = 1;
    while(k < N){
        // std::cout << k << std::endl;
        
        // size reduce
        GridOp<N> G_size_reduce;
        G_size_reduce = GridOp<N>::Identity();
        for(int j = k-1; j >= 0; j--){
            FTYPE mu_kj = mu(k, j);
            FTYPE nu_kj = nu(k, j);

            if(abs_X_aste(j) > abs_Y_aste(j)){
                if(abs(mu_kj) > sqrt2 / 2 || abs(nu_kj) > (1 + sqrt2) / 2){
                    FTYPE x0 = mu_kj - sqrt2 / 2;
                    FTYPE x1 = mu_kj + sqrt2 / 2;
                    FTYPE y0 = nu_kj - (1 + sqrt2) / 2;
                    FTYPE y1 = nu_kj + (1 + sqrt2) / 2;
                    auto solutions = one_dim_grid_problem(x0, x1, y0, y1);
                    ZRoot2 mu_round = solutions[0];
                    ZRoot2 nu_round = conj(mu_round);
                    X.col(k) -= X.col(j) * ZRoot2_to_FTYPE(mu_round);
                    Y.col(k) -= Y.col(j) * ZRoot2_to_FTYPE(nu_round);
                    G_size_reduce(j, k) = -mu_round;
                }
            }else{
                if(abs(mu_kj) > (1 + sqrt2) / 2 || abs(nu_kj) > sqrt2 / 2){
                    FTYPE x0 = mu_kj - (1 + sqrt2) / 2;
                    FTYPE x1 = mu_kj + (1 + sqrt2) / 2;
                    FTYPE y0 = nu_kj - sqrt2 / 2;
                    FTYPE y1 = nu_kj + sqrt2 / 2;
                    auto solutions = one_dim_grid_problem(x0, x1, y0, y1);
                    ZRoot2 mu_round = solutions[0];
                    ZRoot2 nu_round = conj(mu_round);
                    X.col(k) -= X.col(j) * ZRoot2_to_FTYPE(mu_round);
                    Y.col(k) -= Y.col(j) * ZRoot2_to_FTYPE(nu_round);
                    G_size_reduce(j, k) = -mu_round;
                }
            }
            // X_aste = GSO(X);
            // Y_aste = GSO(Y);
        }
        G = G * G_size_reduce;

        FTYPE left = delta * (abs_X_aste(k-1) + abs_Y_aste(k-1)) - mu(k, k-1)*mu(k, k-1) * abs_X_aste(k-1) - nu(k, k-1)*nu(k, k-1) * abs_Y_aste(k-1);
        FTYPE right = abs_X_aste(k) + abs_Y_aste(k);
        if(left <= right){
            k++;
        }else{
            X.col(k-1).swap(X.col(k));
            Y.col(k-1).swap(Y.col(k));
            G.col(k-1).swap(G.col(k));
            X_aste = GSO(X);
            Y_aste = GSO(Y);

            k = std::max(k-1, 1);
        }
    }

    std::cout << "X_reduced = \n" << X << std::endl;
    std::cout << "Y_reduced = \n" << Y << std::endl;

    return G;
}
}