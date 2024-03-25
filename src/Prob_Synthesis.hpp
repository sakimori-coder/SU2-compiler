#pragma once 

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <Eigen/Core>

#include "type.hpp"
#include "quaternion.hpp"
#include "U2_ZOmega.hpp"


namespace SU2_Compiler
{
    /*
    Choi_Jamiolkowski行列を計算
    */
    template <UINT dim>
    Eigen::Matrix<CTYPE, dim*dim, dim*dim> Choi_Jamiolkowski(Eigen::Matrix<CTYPE, dim, dim> U);

    /*
    U(2)のChoi_JamiolkowskiのMagic basis表現を求める
    返り値は実行列
    */
    Eigen::Matrix<FTYPE, 4, 4> CJ_MB(Eigen::Matrix<CTYPE, 2, 2> U);

    /*
    Eigen行列を2次元vectorに変換
    */
    template <typename T, int N, int M>
    std::vector<std::vector<T>> to_2d_vec(const Eigen::Matrix<T, N, M>& A);

    /*
    ユニタリmap Aと確率混合ユニタリmap {B[i], prob[i]}のダイヤモンド距離を求める
    */
    FTYPE distance(quaternion A, std::vector<quaternion> B, std::vector<FTYPE> prob);
    FTYPE distance(quaternion U, std::vector<std::pair<FTYPE, U2_ZOmega>> Prob_Clifford_T);
    FTYPE distance(quaternion U, std::vector<std::pair<FTYPE, std::string>> Prob_Clifford_T);

    /*
    ダイヤモンド距離を最も小さくする確率分布を求める。
    */
    std::pair<FTYPE, std::vector<FTYPE>> get_optimal_prob(std::vector<quaternion> availableU, quaternion targetU);
    
    /*
    与えらたT_countに対して、
    */
    std::pair<std::vector< std::pair<FTYPE, U2_ZOmega> >, std::vector<U2_ZOmega>>
    optimal_prob_unitary(int T_count, quaternion targetU, FTYPE eps, std::vector<U2_ZOmega> pre_availableU_ZOmega = {});

    /*
    
    */
    std::vector<std::pair<FTYPE, std::string>> Prob_Unitary_Synthesis(quaternion targetU, FTYPE eps);
}
