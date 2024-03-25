#pragma once

#include <vector>
#include <utility>
#include <Eigen/Core>


namespace SU2_Compiler
{
    /*
    グラム・シュミットの正規直交化  
    入力 : a_i (i=1,...,M)を正規直交化したいN次元のベクトルとすると、A = (a_1,...,a_M)で与える
    出力 : 正規直交化された基底を並べた行列 V = (v_1,...,v_M)を返す
    */
    template<typename T, int N, int M>
    Eigen::Matrix<T, N, M> Gram_Schmidt(Eigen::Matrix<T, N, M> A);

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
    row_reduction(Eigen::Matrix<T, N, M> A, Eigen::Matrix<T, N, 1> b, std::vector<int> v);

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
    stable_row_reduction(Eigen::Matrix<T, N, M> A, Eigen::Matrix<T, N, 1> b);

    /*
    Ax = bを解く (N <= Mを仮定)

    Returns:
    解空間のM-N本の基底行列(M行(M-N)列)と切片ベクトルのpairを返す
    orthonormalがtrueなら基底は正規直交化して返す。
    */
    template<typename T, int N, int M>
    std::pair<Eigen::Matrix<T, M, M-N>, Eigen::Matrix<T, M, 1>> 
    solve_linear_system(Eigen::Matrix<T, N, M> A, Eigen::Matrix<T, N, 1> b, bool orthonormal=true);
}