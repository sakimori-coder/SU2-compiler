#pragma once

#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "Rings.cpp"
#include "grid_solve.cpp"
#include "search_four_square.cpp"
#include "quaternion.hpp"


template<typename T>
using Pair_ZRoot2 = std::pair<ZRoot2<T>, ZRoot2<T>>;
template<typename T>
using Pair_ZOmega = std::pair<ZOmega<T>, ZOmega<T>>;



template<typename ITYPE, typename FTYPE>
std::vector<quaternion<FTYPE>> enumerate_u_t(quaternion<FTYPE> U, FTYPE eps, int k)
{
    FTYPE a = U.get_a();
    FTYPE b = U.get_b();
    FTYPE c = U.get_c();
    FTYPE d = U.get_d();

    const FTYPE sqrt2 = sqrt((FTYPE)2.0);
    const FTYPE sqrt2k = pow(sqrt2, (FTYPE)k);
    FTYPE y0 =-sqrt2k;
    FTYPE y1 = sqrt2k;


    std::vector<ZRoot2<ITYPE>> X = one_dim_grid_problem<ITYPE, FTYPE>((a-eps)*sqrt2k, (a+eps)*sqrt2k, y0, y1);
    std::vector<ZRoot2<ITYPE>> Y = one_dim_grid_problem<ITYPE, FTYPE>((b-eps)*sqrt2k, (b+eps)*sqrt2k, y0, y1);
    std::vector<ZRoot2<ITYPE>> Z = one_dim_grid_problem<ITYPE, FTYPE>((c-eps)*sqrt2k, (c+eps)*sqrt2k, y0, y1);
    std::vector<ZRoot2<ITYPE>> W = one_dim_grid_problem<ITYPE, FTYPE>((d-eps)*sqrt2k, (d+eps)*sqrt2k, y0, y1);


    FTYPE inv_sqrt2 = 1.0 / sqrt((FTYPE)2.0);
    y0 += inv_sqrt2;
    y1 += inv_sqrt2;
    std::vector<ZRoot2<ITYPE>> X_omega = one_dim_grid_problem<ITYPE, FTYPE>((a-eps)*sqrt2k - inv_sqrt2, (a+eps)*sqrt2k - inv_sqrt2, y0, y1);
    std::vector<ZRoot2<ITYPE>> Y_omega = one_dim_grid_problem<ITYPE, FTYPE>((b-eps)*sqrt2k - inv_sqrt2, (b+eps)*sqrt2k - inv_sqrt2, y0, y1);
    std::vector<ZRoot2<ITYPE>> Z_omega = one_dim_grid_problem<ITYPE, FTYPE>((c-eps)*sqrt2k - inv_sqrt2, (c+eps)*sqrt2k - inv_sqrt2, y0, y1);
    std::vector<ZRoot2<ITYPE>> W_omega = one_dim_grid_problem<ITYPE, FTYPE>((d-eps)*sqrt2k - inv_sqrt2, (d+eps)*sqrt2k - inv_sqrt2, y0, y1);

    // ωの1/√2を消すため, 2(x^2 + y^2 + z^2 + w^2) = 2*2^kを満たすx,y,z,wを求める.
    // 候補点を√2倍したもの 
    std::vector<ZRoot2<ITYPE>> sqrt2_X(X.size());
    std::vector<ZRoot2<ITYPE>> sqrt2_Y(Y.size());
    std::vector<ZRoot2<ITYPE>> sqrt2_Z(Z.size());
    std::vector<ZRoot2<ITYPE>> sqrt2_W(W.size());
    std::vector<ZRoot2<ITYPE>> sqrt2_X_omega(X_omega.size());
    std::vector<ZRoot2<ITYPE>> sqrt2_Y_omega(Y_omega.size());
    std::vector<ZRoot2<ITYPE>> sqrt2_Z_omega(Z_omega.size());
    std::vector<ZRoot2<ITYPE>> sqrt2_W_omega(W_omega.size());

    // √2*(a + b√2) = 2b + a√2
    for(int i = 0; i < X.size(); i++) sqrt2_X[i] = {2*X[i].b, X[i].a};
    for(int i = 0; i < Y.size(); i++) sqrt2_Y[i] = {2*Y[i].b, Y[i].a};
    for(int i = 0; i < Z.size(); i++) sqrt2_Z[i] = {2*Z[i].b, Z[i].a};
    for(int i = 0; i < W.size(); i++) sqrt2_W[i] = {2*W[i].b, W[i].a};
    // √2*(a + b√2 + 1/√2) = (2b + 1) + a√2
    for(int i = 0; i < X_omega.size(); i++) sqrt2_X_omega[i] = {2*X_omega[i].b - 1, X_omega[i].a};
    for(int i = 0; i < Y_omega.size(); i++) sqrt2_Y_omega[i] = {2*Y_omega[i].b - 1, Y_omega[i].a};
    for(int i = 0; i < Z_omega.size(); i++) sqrt2_Z_omega[i] = {2*Z_omega[i].b - 1, Z_omega[i].a};
    for(int i = 0; i < W_omega.size(); i++) sqrt2_W_omega[i] = {2*W_omega[i].b - 1, W_omega[i].a};

    ZRoot2<ITYPE> pow_2_k = (ITYPE)1<<k; 
    ZRoot2<ITYPE> pow_2_k_1 = (ITYPE)1<<(k+1); 
    std::vector<Pair_ZOmega<ITYPE>> solutions_total;   // 全パターンの解の結果 (Z[ω]のpairのvectorであることに注意)

    // u = (a+b√2)+i(c+c√2)をZ[ω]の係数に変換する関数. convert2は√2倍されているものを変換する. 
    auto convert1 = [](ZRoot2<ITYPE> re, ZRoot2<ITYPE> im) -> ZOmega<ITYPE>{ return {-re.b+im.b, im.a, re.b+im.b, re.a}; }; 
    auto convert2 = [](ZRoot2<ITYPE> re, ZRoot2<ITYPE> im) -> ZOmega<ITYPE>{ return {(-re.a+im.a) / 2, im.b, (re.a+im.a) / 2, re.b}; };


    // u = x+iy, t = z+iwと書ける場合(パターン1)
    auto solutions1 = search_four_square_using_vector(X.begin(), X.end(),
                                         Y.begin(), Y.end(),
                                         Z.begin(), Z.end(),
                                         W.begin(), W.end(), pow_2_k);
    for(auto [x, y, z, w] : solutions1) solutions_total.push_back({convert1(*x, *y), convert1(*z, *w)});

    // u = x+iy, t = z+iw+ωと書ける場合(パターン2)
    auto solutions2 = search_four_square_using_vector(sqrt2_X.begin(), sqrt2_X.end(),
                                         sqrt2_Y.begin(), sqrt2_Y.end(),
                                         sqrt2_Z_omega.begin(), sqrt2_Z_omega.end(),
                                         sqrt2_W_omega.begin(), sqrt2_W_omega.end(), pow_2_k_1);
    for(auto [x, y, z, w] : solutions2) solutions_total.push_back({convert2(*x, *y), convert2(*z, *w)});

    // u = x+iy+ω, t = z+iwと書ける場合(パターン3)
    auto solutions3 = search_four_square_using_vector(sqrt2_X_omega.begin(), sqrt2_X_omega.end(),
                                         sqrt2_Y_omega.begin(), sqrt2_Y_omega.end(),
                                         sqrt2_Z.begin(), sqrt2_Z.end(),
                                         sqrt2_W.begin(), sqrt2_W.end(), pow_2_k_1);
    for(auto [x, y, z, w] : solutions3) solutions_total.push_back({convert2(*x, *y), convert2(*z, *w)});

    // u = x+iy+ω, t = z+iw+ωと書ける場合(パターン4)
    auto solutions4 = search_four_square_using_vector(sqrt2_X_omega.begin(), sqrt2_X_omega.end(),
                                         sqrt2_Y_omega.begin(), sqrt2_Y_omega.end(),
                                         sqrt2_Z_omega.begin(), sqrt2_Z_omega.end(),
                                         sqrt2_W_omega.begin(), sqrt2_W_omega.end(), pow_2_k_1);
    for(auto [x, y, z, w] : solutions4) solutions_total.push_back({convert2(*x, *y), convert2(*z, *w)});

    std::vector<quaternion<FTYPE>> ret;
    for(auto [x,y] : solutions_total){
        quaternion<FTYPE> V = to_quaterion<ITYPE, FTYPE>(x, y);
        V = V / V.norm();
        if(distance(U, V) <= eps){
            quaternion<FTYPE> cand1 = U - V;
            quaternion<FTYPE> cand2 = U + V;
            std::cout << std::setprecision(20) << V.norm() << std::endl;
            if(cand1.norm() < cand2.norm()) ret.push_back(V);
            else ret.push_back(-V);
        }
    }

    return ret;
}