#pragma once 

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <Eigen/Core>

#include "type.hpp"
#include "quaternion.hpp"
#include "U2_ZOmega.hpp"


namespace su2compiler
{
    /*
    Choi_Jamiolkowski行列を計算
    */
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
    std::vector< std::pair<FTYPE, U2_ZOmega> >
    optimal_prob_unitary(int T_count, quaternion targetU, FTYPE target_eps);

    /*
    
    */
    std::vector<std::pair<FTYPE, std::string>> Prob_Unitary_Synthesis(quaternion targetU, FTYPE eps);

    int get_T_count(std::vector<std::pair<FTYPE, std::string>> mixed_unitary);

    void set_targetCJUMB();
}
