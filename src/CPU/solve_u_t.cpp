#pragma once

#include <bits/stdc++.h>
#include "Rings.cpp"
#include "grid_solve.cpp"
#include "search_four_square.cpp"


template<typename T>
using Pair_ZRoot2 = std::pair<ZRoot2<T>, ZRoot2<T>>;
template<typename T>
using Pair_ZOmega = std::pair<ZOmega<T>, ZOmega<T>>;



template<typename T>
double distance(
    std::complex<double> u, std::complex<double> t, ZOmega<T> u_approx, ZOmega<T> t_approx, int k)
{
    std::complex<double> diff_u = u - convert(u_approx, k);
    std::complex<double> diff_t = t - convert(t_approx, k);
    double ret1 =  std::sqrt((diff_u * std::conj(diff_u)).real() + (diff_t * std::conj(diff_t)).real());

    diff_u = u + convert(u_approx, k);
    diff_t = t + convert(t_approx, k);
    double ret2 =  std::sqrt((diff_u * std::conj(diff_u)).real() + (diff_t * std::conj(diff_t)).real());

    return std::min(ret1, ret2);
}



template<typename T>
Pair_ZOmega<T> _solve_u_t(std::complex<double> u, std::complex<double> t, double eps, int k)
{
    double a = u.real();
    double b = u.imag();
    double c = t.real();
    double d = t.imag();

    const double sqrt2k = std::pow(std::sqrt(2.0), (double)k);
    double y0 =-sqrt2k;
    double y1 = sqrt2k;

    std::vector<ZRoot2<T>> X = one_dim_grid_problem<T>((a-eps)*sqrt2k, (a+eps)*sqrt2k, y0, y1);
    std::vector<ZRoot2<T>> Y = one_dim_grid_problem<T>((b-eps)*sqrt2k, (b+eps)*sqrt2k, y0, y1);
    std::vector<ZRoot2<T>> Z = one_dim_grid_problem<T>((c-eps)*sqrt2k, (c+eps)*sqrt2k, y0, y1);
    std::vector<ZRoot2<T>> W = one_dim_grid_problem<T>((d-eps)*sqrt2k, (d+eps)*sqrt2k, y0, y1);

    double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    y0 += inv_sqrt2;
    y1 += inv_sqrt2;
    std::vector<ZRoot2<T>> X_omega = one_dim_grid_problem<T>((a-eps)*sqrt2k - inv_sqrt2, (a+eps)*sqrt2k - inv_sqrt2, y0, y1);
    std::vector<ZRoot2<T>> Y_omega = one_dim_grid_problem<T>((b-eps)*sqrt2k - inv_sqrt2, (b+eps)*sqrt2k - inv_sqrt2, y0, y1);
    std::vector<ZRoot2<T>> Z_omega = one_dim_grid_problem<T>((c-eps)*sqrt2k - inv_sqrt2, (c+eps)*sqrt2k - inv_sqrt2, y0, y1);
    std::vector<ZRoot2<T>> W_omega = one_dim_grid_problem<T>((d-eps)*sqrt2k - inv_sqrt2, (d+eps)*sqrt2k - inv_sqrt2, y0, y1);

    // ωの1/√2を消すため, 2(x^2 + y^2 + z^2 + w^2) = 2*2^kを満たすx,y,z,wを求める.
    // 候補点を√2倍したもの 
    std::vector<ZRoot2<T>> sqrt2_X(X.size());
    std::vector<ZRoot2<T>> sqrt2_Y(Y.size());
    std::vector<ZRoot2<T>> sqrt2_Z(Z.size());
    std::vector<ZRoot2<T>> sqrt2_W(W.size());
    std::vector<ZRoot2<T>> sqrt2_X_omega(X_omega.size());
    std::vector<ZRoot2<T>> sqrt2_Y_omega(Y_omega.size());
    std::vector<ZRoot2<T>> sqrt2_Z_omega(Z_omega.size());
    std::vector<ZRoot2<T>> sqrt2_W_omega(W_omega.size());

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

    ZRoot2 pow_2_k = (T)1<<k; 
    ZRoot2 pow_2_k_1 = (T)1<<(k+1); 
    std::vector<Pair_ZOmega<T>> solutions_total;   // 全パターンの解の結果 (Z[ω]のpairのvectorであることに注意)

    // u = (a+b√2)+i(c+c√2)をZ[ω]の係数に変換する関数. convert2は√2倍されているものを変換する. 
    auto convert1 = [](ZRoot2<T> re, ZRoot2<T> im) -> ZOmega<T>{ return {-re.b+im.b, im.a, re.b+im.b, re.a}; }; 
    auto convert2 = [](ZRoot2<T> re, ZRoot2<T> im) -> ZOmega<T>{ return {(-re.a+im.a) / 2, im.b, (re.a+im.a) / 2, re.b}; };


    // u = x+iy, t = z+iwと書ける場合(パターン1)
    auto solutions1 = search_four_square(X.begin(), X.end(),
                                         Y.begin(), Y.end(),
                                         Z.begin(), Z.end(),
                                         W.begin(), W.end(), pow_2_k);
    for(auto [x, y, z, w] : solutions1) solutions_total.push_back({convert1(*x, *y), convert1(*z, *w)});

    // u = x+iy, t = z+iw+ωと書ける場合(パターン2)
    auto solutions2 = search_four_square(sqrt2_X.begin(), sqrt2_X.end(),
                                         sqrt2_Y.begin(), sqrt2_Y.end(),
                                         sqrt2_Z_omega.begin(), sqrt2_Z_omega.end(),
                                         sqrt2_W_omega.begin(), sqrt2_W_omega.end(), pow_2_k_1);
    for(auto [x, y, z, w] : solutions2) solutions_total.push_back({convert2(*x, *y), convert2(*z, *w)});

    // u = x+iy+ω, t = z+iwと書ける場合(パターン3)
    auto solutions3 = search_four_square(sqrt2_X_omega.begin(), sqrt2_X_omega.end(),
                                         sqrt2_Y_omega.begin(), sqrt2_Y_omega.end(),
                                         sqrt2_Z.begin(), sqrt2_Z.end(),
                                         sqrt2_W.begin(), sqrt2_W.end(), pow_2_k_1);
    for(auto [x, y, z, w] : solutions3) solutions_total.push_back({convert2(*x, *y), convert2(*z, *w)});

    // u = x+iy+ω, t = z+iw+ωと書ける場合(パターン4)
    auto solutions4 = search_four_square(sqrt2_X_omega.begin(), sqrt2_X_omega.end(),
                                         sqrt2_Y_omega.begin(), sqrt2_Y_omega.end(),
                                         sqrt2_Z_omega.begin(), sqrt2_Z_omega.end(),
                                         sqrt2_W_omega.begin(), sqrt2_W_omega.end(), pow_2_k_1);
    for(auto [x, y, z, w] : solutions4) solutions_total.push_back({convert2(*x, *y), convert2(*z, *w)});

    double min_error = eps + 1.0;    // 初期値はepsより大きい適当な値
    Pair_ZOmega<T> ret;
    for(auto P : solutions_total){
        ZOmega<T> u_approx = P.first;
        ZOmega<T> t_approx = P.second;
        double error = distance(u, t, u_approx, t_approx, k);
        if(error < min_error){
            min_error = error;
            ret = {u_approx, t_approx};
        }
    }


    if(min_error < eps) return ret;
    else return {};
}


template<typename T>
std::tuple<ZOmega<T>, ZOmega<T>, int> solve_u_t(std::complex<double> u, std::complex<double> t, double eps){
    assert(eps >= 0);
    Pair_ZOmega<T> v = {};
    for(int k = 0; k < 4*std::log2(1/eps); k++){
        // cout << "k " <<  k << endl;
        Pair_ZOmega<T> solution = _solve_u_t<T>(u, t, eps, k);
        if(solution != v) return {solution.first, solution.second, k};
    }
    
    return {};   // warning出るので
}