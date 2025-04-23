#pragma once

#include <iostream>
#include <utility>
#include <vector>

#include "type.hpp"
#include "quaternion.hpp"
#include "linalg.hpp"


namespace SU2_Compiler
{
    /*
    B(u1,ε)とB(u2,ε)が交わるかを判定する。
    出力 : 交わるならtrue, 交わらないならfalse
    */
    bool check_inter_2ball(quaternion U1, quaternion U2, FTYPE eps);

    /*
    B(u1,ε)とB(u2,ε)が交わるある１点を計算する。
    出力 : 交点を１つを出力する
    */
    quaternion cal_inter_2ball(quaternion U1, quaternion U2, FTYPE eps);

    /*
    B(u1,ε)とB(u2,ε)が交わる領域の中で最もtargetUに近い点を計算する。
    */
    quaternion cal_inter_2ball(quaternion U1, quaternion U2, quaternion targetU, FTYPE eps);

    /*
    B(u1,ε)とB(u2,ε)とB(u3,ε)が交わる点を計算する。
    出力 : 交わるならば、交わる2点を出力する。交わらなければ,(0,0,0,0)を返す。
    */
    std::pair<quaternion, quaternion> cal_inter_3ball(quaternion U1, quaternion U2, quaternion U3, FTYPE eps);

    bool check_eps_net(std::vector<quaternion> availableU, quaternion targetU, FTYPE eps);
}