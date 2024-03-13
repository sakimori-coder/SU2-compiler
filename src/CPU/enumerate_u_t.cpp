#pragma once

#include <bits/stdc++.h>
#include <execution>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "Rings.cpp"
#include "grid_solve.cpp"
#include "search_four_square.cpp"
#include "quaternion.hpp"


template<typename T>
using Pair_ZRoot2 = std::pair<ZRoot2<T>, ZRoot2<T>>;
template<typename T>
using Pair_ZOmega = std::pair<ZOmega<T>, ZOmega<T>>;

template <typename T>
T min(T a, T b){
    if(a < b) return a;
    else return b;
}


template <typename T>
std::vector<std::array<T, 4>> subroutine(
    std::vector<T>& X, std::vector<T>& Y, std::vector<T>& Z, std::vector<T>& W,
    std::vector<T>& X_squared, std::vector<T>& Y_squared, std::vector<T>& Z_squared, std::vector<T>& W_squared,
    std::vector<T>& XY, std::vector<T>& ZW, T key)
{
    std::vector<std::array<T, 4>> ret;
    auto solutions = two_points_technique(XY.begin(), XY.end(), ZW.begin(), ZW.end(), key);

    for(auto [xy_squared, zw_squared] : solutions){
        std::set<std::pair<T, T>> xy_ans, zw_ans;
        auto xy_solutions = two_points_technique(X_squared.begin(), X_squared.end(), Y_squared.begin(), Y_squared.end(), xy_squared);
        for(auto [x_square, y_square] : xy_solutions){
            std::vector<T> x_ans, y_ans;
            for(auto x : X) if(x*x == x_square) x_ans.push_back(x);
            for(auto y : Y) if(y*y == y_square) y_ans.push_back(y);
            for(auto x : x_ans) for(auto y : y_ans) xy_ans.insert({x,y});
        }

        auto zw_solutions = two_points_technique(Z_squared.begin(), Z_squared.end(), W_squared.begin(), W_squared.end(), zw_squared);
        for(auto [z_square, w_square]: zw_solutions){
            std::vector<T> z_ans, w_ans;
            for(auto z : Z) if(z*z == z_square) z_ans.push_back(z);
            for(auto w : W) if(w*w == w_square) w_ans.push_back(w);
            for(auto z : z_ans) for(auto w : w_ans) zw_ans.insert({z,w});
        }

        for(auto [x,y] : xy_ans) for(auto [z,w] : zw_ans) ret.push_back({x,y,z,w});
    }

    return ret;
}


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


    // eps *= sqrt2;
    std::vector<ZRoot2<ITYPE>> X = one_dim_grid_problem<ITYPE, FTYPE>((a-eps)*sqrt2k, (a+eps)*sqrt2k, y0, y1);
    std::vector<ZRoot2<ITYPE>> Y = one_dim_grid_problem<ITYPE, FTYPE>((b-eps)*sqrt2k, (b+eps)*sqrt2k, y0, y1);
    std::vector<ZRoot2<ITYPE>> Z = one_dim_grid_problem<ITYPE, FTYPE>((c-eps)*sqrt2k, (c+eps)*sqrt2k, y0, y1);
    std::vector<ZRoot2<ITYPE>> W = one_dim_grid_problem<ITYPE, FTYPE>((d-eps)*sqrt2k, (d+eps)*sqrt2k, y0, y1);

    FTYPE inv_sqrt2 = 1.0 / sqrt((FTYPE)2.0);
    y0 += inv_sqrt2;
    y1 += inv_sqrt2;
    std::vector<ZRoot2<ITYPE>> X_omega = one_dim_grid_problem<ITYPE, FTYPE>((a - eps)*sqrt2k - inv_sqrt2, (a + eps)*sqrt2k - inv_sqrt2, y0, y1);
    std::vector<ZRoot2<ITYPE>> Y_omega = one_dim_grid_problem<ITYPE, FTYPE>((b - eps)*sqrt2k - inv_sqrt2, (b + eps)*sqrt2k - inv_sqrt2, y0, y1);
    std::vector<ZRoot2<ITYPE>> Z_omega = one_dim_grid_problem<ITYPE, FTYPE>((c - eps)*sqrt2k - inv_sqrt2, (c + eps)*sqrt2k - inv_sqrt2, y0, y1);
    std::vector<ZRoot2<ITYPE>> W_omega = one_dim_grid_problem<ITYPE, FTYPE>((d - eps)*sqrt2k - inv_sqrt2, (d + eps)*sqrt2k - inv_sqrt2, y0, y1);
    // eps /= sqrt2;

    const long unsigned int X_size = X.size();
    const long unsigned int Y_size = Y.size();
    const long unsigned int Z_size = Z.size();
    const long unsigned int W_size = W.size();

    const long unsigned int X_omega_size = X_omega.size();
    const long unsigned int Y_omega_size = Y_omega.size();
    const long unsigned int Z_omega_size = Z_omega.size();
    const long unsigned int W_omega_size = W_omega.size();

    // std::cout << X_size << " " << Y_size << " " << Z_size << " " << W_size << std::endl; 
    // std::cout << X_omega_size << " " << Y_omega_size << " " << Z_omega_size << " " << W_omega_size << std::endl; 

    std::vector<ZRoot2<ITYPE>> X_squared(X_size);
    std::vector<ZRoot2<ITYPE>> Y_squared(Y_size);
    std::vector<ZRoot2<ITYPE>> Z_squared(Z_size);
    std::vector<ZRoot2<ITYPE>> W_squared(W_size);
    for(int i = 0; i < X_size; i++) X_squared[i] = X[i] * X[i];
    for(int i = 0; i < Y_size; i++) Y_squared[i] = Y[i] * Y[i];
    for(int i = 0; i < Z_size; i++) Z_squared[i] = Z[i] * Z[i];
    for(int i = 0; i < W_size; i++) W_squared[i] = W[i] * W[i];
    std::sort(std::execution::par, X_squared.begin(), X_squared.end());
    std::sort(std::execution::par, Y_squared.begin(), Y_squared.end());
    std::sort(std::execution::par, Z_squared.begin(), Z_squared.end());
    std::sort(std::execution::par, W_squared.begin(), W_squared.end());


    std::vector<ZRoot2<ITYPE>> X_omega_squared(X_omega_size);
    std::vector<ZRoot2<ITYPE>> Y_omega_squared(Y_omega_size);
    std::vector<ZRoot2<ITYPE>> Z_omega_squared(Z_omega_size);
    std::vector<ZRoot2<ITYPE>> W_omega_squared(W_omega_size);
    // √2倍していないものは使わないので、先に√2倍しておく
    for(auto &x_omega : X_omega) x_omega = {(ITYPE)2*x_omega.b + 1, x_omega.a};
    for(auto &y_omega : Y_omega) y_omega = {(ITYPE)2*y_omega.b + 1, y_omega.a};
    for(auto &z_omega : Z_omega) z_omega = {(ITYPE)2*z_omega.b + 1, z_omega.a};
    for(auto &w_omega : W_omega) w_omega = {(ITYPE)2*w_omega.b + 1, w_omega.a};
    for(int i = 0; i < X_omega_size; i++) X_omega_squared[i] = X_omega[i] * X_omega[i];
    for(int i = 0; i < Y_omega_size; i++) Y_omega_squared[i] = Y_omega[i] * Y_omega[i];
    for(int i = 0; i < Z_omega_size; i++) Z_omega_squared[i] = Z_omega[i] * Z_omega[i];
    for(int i = 0; i < W_omega_size; i++) W_omega_squared[i] = W_omega[i] * W_omega[i];
    std::sort(std::execution::par, X_omega_squared.begin(), X_omega_squared.end());
    std::sort(std::execution::par, Y_omega_squared.begin(), Y_omega_squared.end());
    std::sort(std::execution::par, Z_omega_squared.begin(), Z_omega_squared.end());
    std::sort(std::execution::par, W_omega_squared.begin(), W_omega_squared.end());


    std::vector<ZRoot2<ITYPE>> XY(X_size * Y_size);
    std::vector<ZRoot2<ITYPE>> ZW(Z_size * W_size);
    std::vector<ZRoot2<ITYPE>> XY_omega(X_omega_size * Y_omega_size);
    std::vector<ZRoot2<ITYPE>> ZW_omega(Z_omega_size * W_omega_size);

    for(long long i = 0; i < X_size; i++){
        for(long long j = 0; j < Y_size; j++) XY[i*Y_size + j] = X_squared[i] + Y_squared[j];
    }
    for(long long i = 0; i < Z_size; i++){
        for(long long j = 0; j < W_size; j++) ZW[i*W_size + j] = Z_squared[i] + W_squared[j];
    }
    for(long long i = 0; i < X_omega_size; i++){
        for(long long j = 0; j < Y_omega_size; j++) XY_omega[i*Y_omega_size + j] = X_omega_squared[i] + Y_omega_squared[j];
    }
    for(long long i = 0; i < Z_omega_size; i++){
        for(long long j = 0; j < W_omega_size; j++) ZW_omega[i*W_omega_size + j] = Z_omega_squared[i] + W_omega_squared[j];
    }

#pragma omp parallel sections
    {
#pragma omp section
    std::sort(std::execution::par, XY.begin(), XY.end());
#pragma omp section
    std::sort(std::execution::par, ZW.begin(), ZW.end());
#pragma omp section
    std::sort(std::execution::par, XY_omega.begin(), XY_omega.end());
#pragma omp section
    std::sort(std::execution::par, ZW_omega.begin(), ZW_omega.end());
    }

    ZRoot2<ITYPE> pow_2_k = (ITYPE)1<<k; 
    ZRoot2<ITYPE> pow_2_k_1 = (ITYPE)1<<(k+1); 
    std::vector<Pair_ZOmega<ITYPE>> solutions_total;   // 全パターンの解の結果 (Z[ω]のpairのvectorであることに注意)

    // u = (a+b√2)+i(c+d√2)をZ[ω]の係数に変換する関数. convert2は√2倍されているものを変換する. 
    auto convert1 = [](ZRoot2<ITYPE> re, ZRoot2<ITYPE> im) -> ZOmega<ITYPE>{ return {-re.b+im.b, im.a, re.b+im.b, re.a}; }; 
    auto convert2 = [](ZRoot2<ITYPE> re, ZRoot2<ITYPE> im) -> ZOmega<ITYPE>{ return {(-re.a+im.a) / 2, im.b, (re.a+im.a) / 2, re.b}; };


    auto solutions1 = subroutine(X, Y, Z, W, X_squared, Y_squared, Z_squared, W_squared, XY, ZW, pow_2_k);
    for(auto [x, y, z, w]: solutions1) solutions_total.push_back({convert1(x,y), convert1(z,w)});

    ZRoot2<ITYPE> sqrt2_ZRoot2 = {0,1};
    for(auto &x : X) x = sqrt2_ZRoot2 * x; 
    for(auto &y : Y) y = sqrt2_ZRoot2 * y;
    for(auto &z : Z) z = sqrt2_ZRoot2 * z;
    for(auto &w : W) w = sqrt2_ZRoot2 * w;
    for(auto &x_sq : X_squared) x_sq = x_sq * 2;
    for(auto &y_sq : Y_squared) y_sq = y_sq * 2;
    for(auto &z_sq : Z_squared) z_sq = z_sq * 2;
    for(auto &w_sq : W_squared) w_sq = w_sq * 2;
    for(auto &xy : XY) xy = xy * 2;
    for(auto &zw : ZW) zw = zw * 2;

    auto solutions2 = subroutine(X, Y, Z_omega, W_omega,
                                 X_squared, Y_squared, Z_omega_squared, W_omega_squared, XY, ZW_omega, pow_2_k_1);
    for(auto [x, y, z, w]: solutions2) solutions_total.push_back({convert2(x,y), convert2(z,w)});

    auto solutions3 = subroutine(X_omega, Y_omega, Z, W, 
                                 X_omega_squared, Y_omega_squared, Z_squared, W_squared, 
                                 XY_omega, ZW, pow_2_k_1);
    for(auto [x, y, z, w]: solutions3) solutions_total.push_back({convert2(x,y), convert2(z,w)});

    auto solutions4 = subroutine(X_omega, Y_omega, Z_omega, W_omega,
                                 X_omega_squared, Y_omega_squared, Z_omega_squared, W_omega_squared,
                                XY_omega, ZW_omega, pow_2_k_1);
    for(auto [x, y, z, w]: solutions4) solutions_total.push_back({convert2(x,y), convert2(z,w)});

    std::vector<quaternion<FTYPE>> ret;
    for(auto [x,y] : solutions_total){
        // std::cout << x << " " << y << std::endl;
        quaternion<FTYPE> V = to_quaterion<ITYPE, FTYPE>(x, y);
        V = V / V.norm();

        FTYPE d1 = distance(U,V);
        FTYPE d2 = min(distance_max(U,V), distance_max(U,-V));
        if(d1 < d2/sqrt2) std::cout << "間違ってる" << std::endl;

        // std::cout << d1 << std::endl;

        if(distance(U, V) <= eps){
            quaternion<FTYPE> cand1 = U - V;
            quaternion<FTYPE> cand2 = U + V;

            if(cand1.norm() < cand2.norm()) ret.push_back(V);
            else ret.push_back(-V);
        }
    }

    return ret;
}