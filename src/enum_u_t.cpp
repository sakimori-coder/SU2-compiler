#include "enum_u_t.hpp"

#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <algorithm>
#include <execution>

#include "type.hpp"
#include "rings.hpp"
#include "quaternion.hpp"
#include "U2_ZOmega.hpp"
#include "grid_solver.hpp"
#include "ExactSynthesis.hpp"


#include <chrono>
std::chrono::system_clock::time_point start, end;
double diff_time;
using std::cout;
using std::endl; 



namespace SU2_Compiler
{
    /*
    ソート済み配列A,Bからa_i + b_j = keyとなるi,jを探索する.
    */
    template<class RandomAccessIterator, typename T>
    std::vector<std::pair<T, T>> two_points_technique(
                                RandomAccessIterator first1, RandomAccessIterator last1,
                                RandomAccessIterator first2, RandomAccessIterator last2,
                                const T& key)
    {
        std::vector<std::pair<T, T>> ret;
        static const T zero = 0;
        auto left = first1;
        auto right = last2 - 1;
        for(left = first1; left != last1; left++){
            while(right != first2-1){
                T diff = key - *left - *right;
                if(diff == zero){
                    ret.push_back({*left, *right});
                }
                if(diff < zero) right--;
                else break;
            }
        }

        return ret;
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



    std::pair< std::vector< U2_ZOmega >, std::array<std::vector< ZRoot2 >, 12> >
    enum_u_t(quaternion U, FTYPE eps, int k, int l, const std::array<std::vector< ZRoot2 >, 12>& pre_results)
    {
        quaternion phase(1,0,0,0);
        quaternion sqrt_omega_FTYPE(sqrt_omega.real(), sqrt_omega.imag(), 0, 0);
        for(int i = 0; i < l; i++) phase *= sqrt_omega_FTYPE;
        U = U * phase;
        
        FTYPE a = U.a;
        FTYPE b = U.b;
        FTYPE c = U.c;
        FTYPE d = U.d;

        const FTYPE sqrt2k = pow(sqrt2, (FTYPE)k);
        FTYPE y0 =-sqrt2k;
        FTYPE y1 = sqrt2k;

        FTYPE eps_prime = sqrt(1 - sqrt(1 - eps*eps)) * sqrt2;
        if(isnan(eps_prime)) eps_prime = eps;
        std::vector<ZRoot2> X = one_dim_grid_problem((a-eps_prime)*sqrt2k, (a+eps_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> Y = one_dim_grid_problem((b-eps_prime)*sqrt2k, (b+eps_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> Z = one_dim_grid_problem((c-eps_prime)*sqrt2k, (c+eps_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> W = one_dim_grid_problem((d-eps_prime)*sqrt2k, (d+eps_prime)*sqrt2k, y0, y1);


        y0 += inv_sqrt2;
        y1 += inv_sqrt2;
        std::vector<ZRoot2> X_omega = one_dim_grid_problem((a - eps_prime)*sqrt2k - inv_sqrt2, (a + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> Y_omega = one_dim_grid_problem((b - eps_prime)*sqrt2k - inv_sqrt2, (b + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> Z_omega = one_dim_grid_problem((c - eps_prime)*sqrt2k - inv_sqrt2, (c + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> W_omega = one_dim_grid_problem((d - eps_prime)*sqrt2k - inv_sqrt2, (d + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);


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

        std::vector<ZRoot2> X_squared(X_size);
        std::vector<ZRoot2> Y_squared(Y_size);
        std::vector<ZRoot2> Z_squared(Z_size);
        std::vector<ZRoot2> W_squared(W_size);
        for(int i = 0; i < X_size; i++) X_squared[i] = X[i] * X[i];
        for(int i = 0; i < Y_size; i++) Y_squared[i] = Y[i] * Y[i];
        for(int i = 0; i < Z_size; i++) Z_squared[i] = Z[i] * Z[i];
        for(int i = 0; i < W_size; i++) W_squared[i] = W[i] * W[i];
        std::sort(std::execution::par, X_squared.begin(), X_squared.end());
        std::sort(std::execution::par, Y_squared.begin(), Y_squared.end());
        std::sort(std::execution::par, Z_squared.begin(), Z_squared.end());
        std::sort(std::execution::par, W_squared.begin(), W_squared.end());


        std::vector<ZRoot2> X_omega_squared(X_omega_size);
        std::vector<ZRoot2> Y_omega_squared(Y_omega_size);
        std::vector<ZRoot2> Z_omega_squared(Z_omega_size);
        std::vector<ZRoot2> W_omega_squared(W_omega_size);
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


        std::vector<ZRoot2> diff_X_squared, diff_Y_squared, diff_Z_squared, diff_W_squared,
                                diff_X_omega_squared, diff_Y_omega_squared, diff_Z_omega_squared, diff_W_omega_squared;
        
        std::set_difference(X_squared.begin(), X_squared.end(), pre_results[0].begin(), pre_results[0].end(), std::back_insert_iterator(diff_X_squared));   
        std::set_difference(Y_squared.begin(), Y_squared.end(), pre_results[1].begin(), pre_results[1].end(), std::back_insert_iterator(diff_Y_squared));   
        std::set_difference(Z_squared.begin(), Z_squared.end(), pre_results[2].begin(), pre_results[2].end(), std::back_insert_iterator(diff_Z_squared));   
        std::set_difference(W_squared.begin(), W_squared.end(), pre_results[3].begin(), pre_results[3].end(), std::back_insert_iterator(diff_W_squared));
        std::set_difference(X_omega_squared.begin(), X_omega_squared.end(), pre_results[4].begin(), pre_results[4].end(), std::back_insert_iterator(diff_X_omega_squared));   
        std::set_difference(Y_omega_squared.begin(), Y_omega_squared.end(), pre_results[5].begin(), pre_results[5].end(), std::back_insert_iterator(diff_Y_omega_squared));   
        std::set_difference(Z_omega_squared.begin(), Z_omega_squared.end(), pre_results[6].begin(), pre_results[6].end(), std::back_insert_iterator(diff_Z_omega_squared));   
        std::set_difference(W_omega_squared.begin(), W_omega_squared.end(), pre_results[7].begin(), pre_results[7].end(), std::back_insert_iterator(diff_W_omega_squared));

        const long unsigned int diff_X_size = diff_X_squared.size();
        const long unsigned int diff_Y_size = diff_Y_squared.size();
        const long unsigned int diff_Z_size = diff_Z_squared.size();
        const long unsigned int diff_W_size = diff_W_squared.size();

        const long unsigned int diff_X_omega_size = diff_X_omega_squared.size();
        const long unsigned int diff_Y_omega_size = diff_Y_omega_squared.size();
        const long unsigned int diff_Z_omega_size = diff_Z_omega_squared.size();
        const long unsigned int diff_W_omega_size = diff_W_omega_squared.size();

        // std::cout << X_size << " " << Y_size << " " << Z_size << " " << W_size << std::endl; 
        // std::cout << X_omega_size << " " << Y_omega_size << " " << Z_omega_size << " " << W_omega_size << std::endl;
        // std::cout << diff_X_size << " " << diff_Y_size << " " << diff_Z_size << " " << diff_W_size << std::endl; 
        // std::cout << diff_X_omega_size << " " << diff_Y_omega_size << " " << diff_Z_omega_size << " " << diff_W_omega_size << std::endl;

        std::vector<ZRoot2> diff_XY((diff_X_size * pre_results[1].size()) + (pre_results[0].size() * diff_Y_size) + (diff_X_size * diff_Y_size));
        std::vector<ZRoot2> diff_ZW((diff_Z_size * pre_results[3].size()) + (pre_results[2].size() * diff_W_size) + (diff_Z_size * diff_W_size));
        std::vector<ZRoot2> diff_XY_omega((diff_X_omega_size * pre_results[5].size()) + (pre_results[4].size() * diff_Y_omega_size) + (diff_X_omega_size * diff_Y_omega_size));
        std::vector<ZRoot2> diff_ZW_omega((diff_Z_omega_size * pre_results[7].size()) + (pre_results[6].size() * diff_W_omega_size) + (diff_Z_omega_size * diff_W_omega_size));


    #pragma omp parallel sections
        {
    #pragma omp section
            {
                auto diff_XY_ite = diff_XY.begin();
                for(auto diff_x_sq : diff_X_squared){
                    for(auto y_sq : pre_results[1]) {*diff_XY_ite = diff_x_sq + y_sq; diff_XY_ite++;}
                }
                for(auto x_sq : pre_results[0]){
                    for(auto diff_y_sq : diff_Y_squared) {*diff_XY_ite = x_sq + diff_y_sq; diff_XY_ite++;}
                }
                for(auto diff_x_sq : diff_X_squared){
                    for(auto diff_y_sq : diff_Y_squared) {*diff_XY_ite = diff_x_sq + diff_y_sq; diff_XY_ite++;}
                }
            }

    #pragma omp section
            {            
                auto diff_ZW_ite = diff_ZW.begin();
                for(auto diff_z_sq : diff_Z_squared){
                    for(auto w_sq : pre_results[3]) {*diff_ZW_ite = diff_z_sq + w_sq; diff_ZW_ite++;}
                }
                for(auto z_sq : pre_results[2]){
                    for(auto diff_w_sq : diff_W_squared) {*diff_ZW_ite = z_sq + diff_w_sq; diff_ZW_ite++;}
                }
                for(auto diff_z_sq : diff_Z_squared){
                    for(auto diff_w_sq : diff_W_squared) {*diff_ZW_ite = diff_z_sq + diff_w_sq; ++diff_ZW_ite;}
                }
            }

    #pragma omp section
            {
                auto diff_XY_omega_ite = diff_XY_omega.begin();
                for(auto diff_x_omega_sq : diff_X_omega_squared){
                    for(auto y_omega_sq : pre_results[5]) {*diff_XY_omega_ite = diff_x_omega_sq + y_omega_sq; diff_XY_omega_ite++;}
                }
                for(auto x_omega_sq : pre_results[4]){
                    for(auto diff_y_omega_sq : diff_Y_omega_squared) {*diff_XY_omega_ite = x_omega_sq + diff_y_omega_sq; diff_XY_omega_ite++;}
                }
                for(auto diff_x_omega_sq : diff_X_omega_squared){
                    for(auto diff_y_omega_sq : diff_Y_omega_squared) {*diff_XY_omega_ite = diff_x_omega_sq + diff_y_omega_sq; diff_XY_omega_ite++;}
                }
            }

    #pragma omp section
            {
                auto diff_ZW_omega_ite = diff_ZW_omega.begin();
                for(auto diff_z_omega_sq : diff_Z_omega_squared){
                    for(auto w_omega_sq : pre_results[7]) {*diff_ZW_omega_ite = diff_z_omega_sq + w_omega_sq; diff_ZW_omega_ite++;}
                }
                for(auto z_omega_sq : pre_results[6]){
                    for(auto diff_w_omega_sq : diff_W_omega_squared) {*diff_ZW_omega_ite = z_omega_sq + diff_w_omega_sq; diff_ZW_omega_ite++;}
                }
                for(auto diff_z_omega_sq : diff_Z_omega_squared){
                    for(auto diff_w_omega_sq : diff_W_omega_squared) {*diff_ZW_omega_ite = diff_z_omega_sq + diff_w_omega_sq; diff_ZW_omega_ite++;}
                }
            }
        }

    #pragma omp parallel sections
        {
    #pragma omp section
            std::sort(std::execution::par, diff_XY.begin(), diff_XY.end());
    #pragma omp section
            std::sort(std::execution::par, diff_ZW.begin(), diff_ZW.end());
    #pragma omp section
            std::sort(std::execution::par, diff_XY_omega.begin(), diff_XY_omega.end());
    #pragma omp section
            std::sort(std::execution::par, diff_ZW_omega.begin(), diff_ZW_omega.end());
        }

        std::vector<ZRoot2> XY(pre_results[8].size() + diff_XY.size());
        std::vector<ZRoot2> ZW(pre_results[9].size() + diff_ZW.size());
        std::vector<ZRoot2> XY_omega(pre_results[10].size() + diff_XY_omega.size());
        std::vector<ZRoot2> ZW_omega(pre_results[11].size() + diff_ZW_omega.size());


    #pragma omp parallel sections
        {
    #pragma omp section
            std::merge(diff_XY.begin(), diff_XY.end(), pre_results[8].begin(), pre_results[8].end(), XY.begin());
    #pragma omp section
            std::merge(diff_ZW.begin(), diff_ZW.end(), pre_results[9].begin(), pre_results[9].end(), ZW.begin());
    #pragma omp section
            std::merge(diff_XY_omega.begin(), diff_XY_omega.end(), pre_results[10].begin(), pre_results[10].end(), XY_omega.begin());
    #pragma omp section
            std::merge(diff_ZW_omega.begin(), diff_ZW_omega.end(), pre_results[11].begin(), pre_results[11].end(), ZW_omega.begin());
        }

        ZRoot2 pow_2_k = (ITYPE)1<<k; 
        ZRoot2 pow_2_k_1 = (ITYPE)1<<(k+1); 
        std::vector<U2_ZOmega> solutions_total;   // 全パターンの解の結果 (Z[ω]のpairのvectorであることに注意)

        // u = (a+b√2)+i(c+d√2)をZ[ω]の係数に変換する関数. convert2は√2倍されているものを変換する. 
        auto convert1 = [](ZRoot2 re, ZRoot2 im) -> ZOmega{ return {-re.b+im.b, im.a, re.b+im.b, re.a}; }; 
        auto convert2 = [](ZRoot2 re, ZRoot2 im) -> ZOmega{ return {(-re.a+im.a) / 2, im.b, (re.a+im.a) / 2, re.b}; };


        auto XYZW1 = subroutine(X, Y, Z, W, X_squared, Y_squared, Z_squared, W_squared, XY, ZW, pow_2_k);
        for(auto [x, y, z, w]: XYZW1) solutions_total.push_back({convert1(x,y), convert1(z,w), l, k});

        ZRoot2 sqrt2_ZRoot2 = {0,1};
        for(auto &x : X) x = sqrt2_ZRoot2 * x; 
        for(auto &y : Y) y = sqrt2_ZRoot2 * y;
        for(auto &z : Z) z = sqrt2_ZRoot2 * z;
        for(auto &w : W) w = sqrt2_ZRoot2 * w;
        for(auto &x_sq : X_squared) x_sq = x_sq * 2;
        for(auto &y_sq : Y_squared) y_sq = y_sq * 2;
        for(auto &z_sq : Z_squared) z_sq = z_sq * 2;
        for(auto &w_sq : W_squared) w_sq = w_sq * 2;

    #pragma omp parallel for
        for(auto &xy : XY) xy = xy * 2;
    #pragma omp parallel for
        for(auto &zw : ZW) zw = zw * 2;

        std::vector<U2_ZOmega> solutions2, solutions3, solutions4;
    #pragma omp parallel sections
        {
    #pragma omp section
            {
                auto XYZW2 = subroutine(X, Y, Z_omega, W_omega,
                                        X_squared, Y_squared, Z_omega_squared, W_omega_squared, XY, ZW_omega, pow_2_k_1);
                for(auto [x, y, z, w]: XYZW2) solutions2.push_back({convert2(x,y), convert2(z,w), l, k});
            }
    #pragma omp section
            {
                auto XYZW3 = subroutine(X_omega, Y_omega, Z, W, 
                                        X_omega_squared, Y_omega_squared, Z_squared, W_squared, 
                                        XY_omega, ZW, pow_2_k_1);
                for(auto [x, y, z, w]: XYZW3) solutions3.push_back({convert2(x,y), convert2(z,w), l, k});
            }
    #pragma omp section
            {
                auto XYZW4 = subroutine(X_omega, Y_omega, Z_omega, W_omega,
                                        X_omega_squared, Y_omega_squared, Z_omega_squared, W_omega_squared,
                                        XY_omega, ZW_omega, pow_2_k_1);
                for(auto [x, y, z, w]: XYZW4) solutions4.push_back({convert2(x,y), convert2(z,w), l, k});
            }
        }

        for(auto ele : solutions2) solutions_total.push_back(ele);
        for(auto ele : solutions3) solutions_total.push_back(ele);
        for(auto ele : solutions4) solutions_total.push_back(ele);


        for(auto &x_sq : X_squared) x_sq = {x_sq.a / 2, x_sq.b / 2};
        for(auto &y_sq : Y_squared) y_sq = {y_sq.a / 2, y_sq.b / 2};
        for(auto &z_sq : Z_squared) z_sq = {z_sq.a / 2, z_sq.b / 2};
        for(auto &w_sq : W_squared) w_sq = {w_sq.a / 2, w_sq.b / 2};

    #pragma omp parallel for
        for(auto &xy : XY) xy = {xy.a / 2, xy.b / 2};
    #pragma omp parallel for
        for(auto &zw : ZW) zw = {zw.a / 2, zw.b / 2};

        std::vector<U2_ZOmega> ret;
        for(auto ele : solutions_total){
            quaternion V = to_quaternion(ele);

            // std::cout << U_copy << std::endl;
            // std::cout << V << std::endl;

            FTYPE d1 = distance(U, V);
            // std::cout << d1 << std::endl;
            // FTYPE d2 = min(distance_max(U,V), distance_max(U,-V));
            // if(d1 < d2/sqrt2) std::cout << "間違ってる" << std::endl;

            // std::cout << d1 << std::endl;

            if(distance(U, V) <= eps){
                // std::string str = Exact_synthesis(u, t, k);
                // cout << str << endl;
                // cout << std::count(str.begin(), str.end(), 'T') << endl;

                // quaternion<FTYPE> cand1 = U - V;
                // quaternion<FTYPE> cand2 = U + V;

                // if(cand1.norm() < cand2.norm()) ret.push_back({u, t});
                // else ret.push_back({-u, -t});

                bool flag = true;
                for(auto &U : ret) if(distance(to_quaternion(U), V) < 1e-20) flag = false;   // 同じユニタリは追加しない
                
                if(flag)ret.push_back(ele);
            }
        }

        auto comp = [](const U2_ZOmega& l, const U2_ZOmega& r){
            return l.u < r.u || ((l.u == r.u) && l.t < r.t);
        };
        std::sort(ret.begin(), ret.end(), comp);
        ret.erase(std::unique(ret.begin(), ret.end()), ret.end());

        return {ret, {X_squared, Y_squared, Z_squared, W_squared,
                    X_omega_squared, Y_omega_squared, Z_omega_squared, W_omega_squared, 
                    XY, ZW, XY_omega, ZW_omega}};
    }
}