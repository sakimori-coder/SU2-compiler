#pragma once

#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <Eigen/Core>


/*
グラム・シュミットの正規直交化  
入力 : a_i (i=1,...,M)を正規直交化したいN次元のベクトルとすると、A = (a_1,...,a_M)で与える
出力 : 正規直交化された基底を並べた行列 V = (v_1,...,v_M)を返す
*/
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



/*
Ax=bの掃き出し法　計算量O(NM*min(N,M))

Parameters:
A : N行M列の行列
b : N次元のベクトル
v : 簡約化する列を指定するvector。サイズはmin(N, M)である必要がある。

Return:
簡約化されたAとbを返す。

Example
Aが4行2列、bが4次元ベクトルであるとする。v=(0,3)である場合は、
[[1 a 0 c],    [[b1],
 [0 b 1 d]]     [b2]]    のような簡約化された行列とベクトルが返される。
*/ 
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



/*
数値的な安定性を考慮したAx=bの掃き出し法　計算量O(NM*min(N,M))

Parameters:
A : N行M列の行列
b : N次元のベクトル

Return:
簡約化されたAとbを返す。また、簡約化した列を表すvectorを返す。

Example
Aが4行2列、bが4次元ベクトルであるとする。selected_rows={0,3}である場合は、
[[1 a c 0],    [[b1],
 [0 b d 1]]     [b2]]    のような簡約化された行列とベクトルが返される。
*/ 
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

/*
Ax = bを解く (N <= Mを仮定)

Returns:
解空間のM-N本の基底行列(M行(M-N)列)と切片ベクトルのpairを返す
orthonormalがtrueなら基底は正規直交化して返す。
*/
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