#pragma once

#include <iostream>
#include <vector>
#include <utility>
#include <numeric>
#include <Eigen/Core>


namespace SU2_Compiler
{
    template<typename T, int N, int M>
    Eigen::Matrix<T, N, M> Gram_Schmidt(Eigen::Matrix<T, N, M> A)
    {
        Eigen::Matrix<T, N, M> V;
        for(int i = 0; i < M; i++){
            Eigen::Matrix<T, N, 1> v = A.col(i);
            for(int j = 0; j < i; j++) v -= V.col(j).dot(A.col(i)) * V.col(j);
            v = v / v.norm();
            V.col(i) = v;
        }

        return V;
    }


    template<typename T, int N, int M>
    std::pair<Eigen::Matrix<T, N, M>, Eigen::Matrix<T, N, 1>>
    row_reduction(Eigen::Matrix<T, N, M> A, Eigen::Matrix<T, N, 1> b, std::vector<int> v)
    {
        for(int k = 0; k < std::min(N, M); k++){
            T a_kvk = A(k,v[k]);
            for(int j = 0; j < M; j++) A(k,j) /= a_kvk;
            b(k) /= a_kvk;

            for(int i = 0; i < N; i++){
                if(i == k) continue;
                T a_ivk = A(i,v[k]);
                for(int j = 0; j < M; j++) A(i,j) -= a_ivk*A(k, j);
                b(i) -= a_ivk*b(k);
            }
        }

        return {A, b};
    }


    template<typename T, int N, int M>
    std::tuple<Eigen::Matrix<T, N, M>, Eigen::Matrix<T, N, 1>, std::vector<int>>
    stable_row_reduction(Eigen::Matrix<T, N, M> A, Eigen::Matrix<T, N, 1> b)
    {
        std::vector<int> selected_rows;   // 簡約化した列を記録
        for(int k = 0; k < std::min(N, M); k++){
            std::vector<T> v(M);  
            for(int i = 0; i < M; i++) v[i] = abs(A(k,i));  // k行目の絶対値を格納する 
            int l = std::distance(v.begin(), std::max_element(v.begin(), v.end()));   // k行目の絶対値が最も大きい列を簡約化する。
                                                                                    // 既に簡約化された列は0が格納されているはずなので選択されない。 
            selected_rows.push_back(l);

            T a_kl = A(k,l);
            for(int j = 0; j < M; j++) A(k,j) /= a_kl;
            b(k) /= a_kl;

            for(int i = 0; i < N; i++){
                if(i == k) continue;
                T a_il = A(i,l);
                for(int j = 0; j < M; j++) A(i,j) -= a_il*A(k, j);
                b(i) -= a_il*b(k);
            }
        }

        return {A, b, selected_rows};
    }


    template<typename T, int N, int M>
    std::pair<Eigen::Matrix<T, M, M-N>, Eigen::Matrix<T, M, 1>> solve_linear_system(Eigen::Matrix<T, N, M> A, Eigen::Matrix<T, N, 1> b, bool orthonormal=true){

        auto[Ar, br, selected_rows] = stable_row_reduction<T, N, M>(A, b);

        std::vector<int> copy_selected_rows = selected_rows;
        std::sort(copy_selected_rows.begin(), copy_selected_rows.end());
        std::vector<int> all_rows(M);
        std::iota(all_rows.begin(), all_rows.end(), 0);   // all_rows = {0, ... , M-1}
        std::vector<int> not_selected_rows;   // 簡約化されていない列の添字集合
        std::set_difference(all_rows.begin(), all_rows.end(), copy_selected_rows.begin(), copy_selected_rows.end(), std::back_inserter(not_selected_rows));

        Eigen::Matrix<T, M, M-N> V;
        for(int i = 0; i < M-N; i++){
            Eigen::Matrix<T, M, 1> v;   // 解空間の基底
            v = Eigen::Matrix<T, M, 1>::Zero(M);

            for(int j = 0; j < N; j++) v[selected_rows[j]] = -Ar(j, not_selected_rows[i]);
            v[not_selected_rows[i]] = 1;

            V.col(i) = v;
        }

        Eigen::Matrix<T, M, 1> inter;
        inter = Eigen::Matrix<T, M, 1>::Zero(M);
        for(int j = 0; j < N; j++) inter[selected_rows[j]] = br(j);

        if(orthonormal) V = Gram_Schmidt<T, M, M-N>(V);

        return {V, inter};
    }
}