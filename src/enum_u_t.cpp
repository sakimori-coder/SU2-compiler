#include "enum_u_t.hpp"

#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <execution>
#include <chrono>
#include <tbb/concurrent_hash_map.h>

#include "type.hpp"
#include "rings.hpp"
#include "quaternion.hpp"
#include "U2_ZOmega.hpp"
#include "grid_solver.hpp"
#include "ExactSynthesis.hpp"
#include "HashTable.hpp"


#include <chrono>
using std::cout;
using std::endl; 



namespace SU2_Compiler
{
    std::chrono::system_clock::time_point start, end;
    double time;

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
        start = std::chrono::system_clock::now();

        quaternion phase(1,0,0,0);
        quaternion sqrt_omega_FTYPE(sqrt_omega.real(), sqrt_omega.imag(), 0, 0);
        for(int i = 0; i < l; i++) phase *= sqrt_omega_FTYPE;
        quaternion U_prime = U * phase;
        
        FTYPE a = U_prime.a;
        FTYPE b = U_prime.b;
        FTYPE c = U_prime.c;
        FTYPE d = U_prime.d;

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
        
        std::set_difference(X_squared.begin(), X_squared.end(), pre_results[0].begin(), pre_results[0].end(), std::back_inserter(diff_X_squared));   
        std::set_difference(Y_squared.begin(), Y_squared.end(), pre_results[1].begin(), pre_results[1].end(), std::back_inserter(diff_Y_squared));   
        std::set_difference(Z_squared.begin(), Z_squared.end(), pre_results[2].begin(), pre_results[2].end(), std::back_inserter(diff_Z_squared));   
        std::set_difference(W_squared.begin(), W_squared.end(), pre_results[3].begin(), pre_results[3].end(), std::back_inserter(diff_W_squared));
        std::set_difference(X_omega_squared.begin(), X_omega_squared.end(), pre_results[4].begin(), pre_results[4].end(), std::back_inserter(diff_X_omega_squared));   
        std::set_difference(Y_omega_squared.begin(), Y_omega_squared.end(), pre_results[5].begin(), pre_results[5].end(), std::back_inserter(diff_Y_omega_squared));   
        std::set_difference(Z_omega_squared.begin(), Z_omega_squared.end(), pre_results[6].begin(), pre_results[6].end(), std::back_inserter(diff_Z_omega_squared));   
        std::set_difference(W_omega_squared.begin(), W_omega_squared.end(), pre_results[7].begin(), pre_results[7].end(), std::back_inserter(diff_W_omega_squared));

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

        cout << diff_XY.size() << endl;
        cout << diff_ZW.size() << endl;
        cout << diff_XY_omega.size() << endl;
        cout << diff_ZW_omega.size() << endl;
        start = std::chrono::system_clock::now();


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

        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        cout << "時間1 " << time << endl;

        start = std::chrono::system_clock::now();

        std::vector<ZRoot2> XY(pre_results[8].size() + diff_XY.size());
        std::vector<ZRoot2> ZW(pre_results[9].size() + diff_ZW.size());
        std::vector<ZRoot2> XY_omega(pre_results[10].size() + diff_XY_omega.size());
        std::vector<ZRoot2> ZW_omega(pre_results[11].size() + diff_ZW_omega.size());

    #pragma omp parallel sections
        {
    #pragma omp section
            std::merge(std::execution::par, diff_XY.begin(), diff_XY.end(), pre_results[8].begin(), pre_results[8].end(), XY.begin());
    #pragma omp section
            std::merge(std::execution::par, diff_ZW.begin(), diff_ZW.end(), pre_results[9].begin(), pre_results[9].end(), ZW.begin());
    #pragma omp section
            std::merge(std::execution::par, diff_XY_omega.begin(), diff_XY_omega.end(), pre_results[10].begin(), pre_results[10].end(), XY_omega.begin());
    #pragma omp section
            std::merge(std::execution::par, diff_ZW_omega.begin(), diff_ZW_omega.end(), pre_results[11].begin(), pre_results[11].end(), ZW_omega.begin());
        }

        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        cout << "時間2 " << time << endl;

        ZRoot2 pow_2_k = (ITYPE)1<<k; 
        ZRoot2 pow_2_k_1 = (ITYPE)1<<(k+1); 
        std::vector<U2_ZOmega> solutions_total;   // 全パターンの解の結果 (Z[ω]のpairのvectorであることに注意)

        // u = (a+b√2)+i(c+d√2)をZ[ω]の係数に変換する関数. convert2は√2倍されているものを変換する. 
        auto convert1 = [](ZRoot2 re, ZRoot2 im) -> ZOmega{ return {-re.b+im.b, im.a, re.b+im.b, re.a}; }; 
        auto convert2 = [](ZRoot2 re, ZRoot2 im) -> ZOmega{ return {(-re.a+im.a) / 2, im.b, (re.a+im.a) / 2, re.b}; };

        
        start = start = std::chrono::system_clock::now();

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

        cout << "候補数 " << solutions_total.size() << endl;

        std::vector<U2_ZOmega> ret;
        for(auto ele : solutions_total){
            quaternion V = to_quaternion(ele);
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

        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        cout << "時間3 " << time << endl;

        return {ret, {X_squared, Y_squared, Z_squared, W_squared,
                    X_omega_squared, Y_omega_squared, Z_omega_squared, W_omega_squared, 
                    XY, ZW, XY_omega, ZW_omega}};
    }


    std::vector<U2_ZOmega>
    enum_u_t
    (quaternion U, FTYPE eps, FTYPE eps_pre, int k, int l, std::vector<ZRoot2>& XY, std::vector<ZRoot2>& ZW, std::vector<U2_ZOmega>& candidates)
    {
        quaternion phase(1,0,0,0);
        quaternion sqrt_omega_FTYPE(sqrt_omega.real(), sqrt_omega.imag(), 0, 0);
        for(int i = 0; i < l; i++) phase *= sqrt_omega_FTYPE;
        quaternion U_prime = U * phase;
        
        
        FTYPE a = U_prime.a;
        FTYPE b = U_prime.b;
        FTYPE c = U_prime.c;
        FTYPE d = U_prime.d;

        const FTYPE sqrt2k = pow(sqrt2, (FTYPE)k);
        FTYPE y0 =-sqrt2k;
        FTYPE y1 = sqrt2k;

        FTYPE eps_prime = sqrt(1 - sqrt(1 - eps * eps)) * sqrt2;
        if(isnan(eps_prime)) eps_prime = eps;

        FTYPE eps_pre_prime = sqrt(1 - sqrt(1 - eps_pre * eps_pre)) * sqrt2;
        if(isnan(eps_pre_prime)) eps_pre_prime = eps_pre;

        std::vector<ZRoot2> X_total = one_dim_grid_problem((a - eps_prime)*sqrt2k, (a + eps_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> Y_total = one_dim_grid_problem((b - eps_prime)*sqrt2k, (b + eps_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> Z_total = one_dim_grid_problem((c - eps_prime)*sqrt2k, (c + eps_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> W_total = one_dim_grid_problem((d - eps_prime)*sqrt2k, (d + eps_prime)*sqrt2k, y0, y1);

        std::vector<ZRoot2> X_pre = one_dim_grid_problem((a - eps_pre_prime)*sqrt2k, (a + eps_pre_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> Y_pre = one_dim_grid_problem((b - eps_pre_prime)*sqrt2k, (b + eps_pre_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> Z_pre = one_dim_grid_problem((c - eps_pre_prime)*sqrt2k, (c + eps_pre_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> W_pre = one_dim_grid_problem((d - eps_pre_prime)*sqrt2k, (d + eps_pre_prime)*sqrt2k, y0, y1);

        // ((√2x)^2 + (√2x)^2 + (√2x)^2 + (√2x)^2) = 2^{k+1}を満たす(x,y,z,w)を求めるので、√2倍する
        const ZRoot2 sqrt2_ZRoot2 = {0,1};
        for(auto &x : X_total) x *= sqrt2_ZRoot2;
        for(auto &y : Y_total) y *= sqrt2_ZRoot2;
        for(auto &z : Z_total) z *= sqrt2_ZRoot2;
        for(auto &w : W_total) w *= sqrt2_ZRoot2;
        std::sort(std::execution::par, X_total.begin(), X_total.end());
        std::sort(std::execution::par, Y_total.begin(), Y_total.end());
        std::sort(std::execution::par, Z_total.begin(), Z_total.end());
        std::sort(std::execution::par, W_total.begin(), W_total.end());

        for(auto &x : X_pre) x *= sqrt2_ZRoot2;
        for(auto &y : Y_pre) y *= sqrt2_ZRoot2;
        for(auto &z : Z_pre) z *= sqrt2_ZRoot2;
        for(auto &w : W_pre) w *= sqrt2_ZRoot2;
        std::sort(std::execution::par, X_pre.begin(), X_pre.end());
        std::sort(std::execution::par, Y_pre.begin(), Y_pre.end());
        std::sort(std::execution::par, Z_pre.begin(), Z_pre.end());
        std::sort(std::execution::par, W_pre.begin(), W_pre.end());

        std::vector<ZRoot2> X_new, Y_new, Z_new, W_new;
        std::set_difference(std::execution::par, X_total.begin(), X_total.end(), X_pre.begin(), X_pre.end(), std::back_inserter(X_new));
        std::set_difference(std::execution::par, Y_total.begin(), Y_total.end(), Y_pre.begin(), Y_pre.end(), std::back_inserter(Y_new));
        std::set_difference(std::execution::par, Z_total.begin(), Z_total.end(), Z_pre.begin(), Z_pre.end(), std::back_inserter(Z_new));
        std::set_difference(std::execution::par, W_total.begin(), W_total.end(), W_pre.begin(), W_pre.end(), std::back_inserter(W_new));

        y0 += inv_sqrt2;
        y1 += inv_sqrt2;
        std::vector<ZRoot2> X_omega_total = one_dim_grid_problem((a - eps_prime)*sqrt2k - inv_sqrt2, (a + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> Y_omega_total = one_dim_grid_problem((b - eps_prime)*sqrt2k - inv_sqrt2, (b + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> Z_omega_total = one_dim_grid_problem((c - eps_prime)*sqrt2k - inv_sqrt2, (c + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> W_omega_total = one_dim_grid_problem((d - eps_prime)*sqrt2k - inv_sqrt2, (d + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);

        std::vector<ZRoot2> X_omega_pre = one_dim_grid_problem((a - eps_pre_prime)*sqrt2k - inv_sqrt2, (a + eps_pre_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> Y_omega_pre = one_dim_grid_problem((b - eps_pre_prime)*sqrt2k - inv_sqrt2, (b + eps_pre_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> Z_omega_pre = one_dim_grid_problem((c - eps_pre_prime)*sqrt2k - inv_sqrt2, (c + eps_pre_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> W_omega_pre = one_dim_grid_problem((d - eps_pre_prime)*sqrt2k - inv_sqrt2, (d + eps_pre_prime)*sqrt2k - inv_sqrt2, y0, y1);

        for(auto &x : X_omega_total) x = (x * sqrt2_ZRoot2) + 1;
        for(auto &y : Y_omega_total) y = (y * sqrt2_ZRoot2) + 1;
        for(auto &z : Z_omega_total) z = (z * sqrt2_ZRoot2) + 1;
        for(auto &w : W_omega_total) w = (w * sqrt2_ZRoot2) + 1;
        std::sort(std::execution::par, X_omega_total.begin(), X_omega_total.end());
        std::sort(std::execution::par, Y_omega_total.begin(), Y_omega_total.end());
        std::sort(std::execution::par, Z_omega_total.begin(), Z_omega_total.end());
        std::sort(std::execution::par, W_omega_total.begin(), W_omega_total.end());

        for(auto &x : X_omega_pre) x = (x * sqrt2_ZRoot2) + 1;
        for(auto &y : Y_omega_pre) y = (y * sqrt2_ZRoot2) + 1;
        for(auto &z : Z_omega_pre) z = (z * sqrt2_ZRoot2) + 1;
        for(auto &w : W_omega_pre) w = (w * sqrt2_ZRoot2) + 1;
        std::sort(std::execution::par, X_omega_pre.begin(), X_omega_pre.end());
        std::sort(std::execution::par, Y_omega_pre.begin(), Y_omega_pre.end());
        std::sort(std::execution::par, Z_omega_pre.begin(), Z_omega_pre.end());
        std::sort(std::execution::par, W_omega_pre.begin(), W_omega_pre.end());

        std::vector<ZRoot2> X_omega_new, Y_omega_new, Z_omega_new, W_omega_new;
        std::set_difference(std::execution::par, X_omega_total.begin(), X_omega_total.end(), X_omega_pre.begin(), X_omega_pre.end(), std::back_inserter(X_omega_new));
        std::set_difference(std::execution::par, Y_omega_total.begin(), Y_omega_total.end(), Y_omega_pre.begin(), Y_omega_pre.end(), std::back_inserter(Y_omega_new));
        std::set_difference(std::execution::par, Z_omega_total.begin(), Z_omega_total.end(), Z_omega_pre.begin(), Z_omega_pre.end(), std::back_inserter(Z_omega_new));
        std::set_difference(std::execution::par, W_omega_total.begin(), W_omega_total.end(), W_omega_pre.begin(), W_omega_pre.end(), std::back_inserter(W_omega_new));

        const long unsigned int X_new_size = X_new.size();
        const long unsigned int Y_new_size = Y_new.size();
        const long unsigned int Z_new_size = Z_new.size();
        const long unsigned int W_new_size = W_new.size();

        const long unsigned int X_pre_size = X_pre.size();
        const long unsigned int Y_pre_size = Y_pre.size();
        const long unsigned int Z_pre_size = Z_pre.size();
        const long unsigned int W_pre_size = W_pre.size();

        const long unsigned int X_omega_new_size = X_omega_new.size();
        const long unsigned int Y_omega_new_size = Y_omega_new.size();
        const long unsigned int Z_omega_new_size = Z_omega_new.size();
        const long unsigned int W_omega_new_size = W_omega_new.size();

        const long unsigned int X_omega_pre_size = X_omega_pre.size();
        const long unsigned int Y_omega_pre_size = Y_omega_pre.size();
        const long unsigned int Z_omega_pre_size = Z_omega_pre.size();
        const long unsigned int W_omega_pre_size = W_omega_pre.size();

        // std::cout << X_new_size + X_pre_size << " " << Y_new_size + Y_pre_size << " " << Z_new_size + Z_pre_size << " " << W_new_size + W_pre_size << std::endl; 
        // std::cout << X_omega_new_size << " " << Y_omega_new_size << " " << Z_omega_new_size << " " << W_omega_new_size << std::endl; 
        // std::cout << X_omega_pre_size << " " << Y_omega_pre_size << " " << Z_omega_pre_size << " " << W_omega_pre_size << std::endl; 
        // std::cout << X_omega_new_size + X_omega_pre_size << " " << Y_omega_new_size + Y_omega_pre_size << " " << Z_omega_new_size + Z_omega_pre_size << " " << W_omega_new_size + W_omega_pre_size << std::endl; 

        // 事前にソートしておくと後のソートの交換回数が減り高速になるかも（ならなかった）
        // auto comp = [](const ZRoot2& x, const ZRoot2& y){return x*x < y*y;};
        // std::sort(std::execution::par, X_new.begin(), X_new.end(), comp);
        // std::sort(std::execution::par, Y_new.begin(), Y_new.end(), comp);
        // std::sort(std::execution::par, Z_new.begin(), Z_new.end(), comp);
        // std::sort(std::execution::par, W_new.begin(), W_new.end(), comp);
        // std::sort(std::execution::par, X_pre.begin(), X_pre.end(), comp);
        // std::sort(std::execution::par, Y_pre.begin(), Y_pre.end(), comp);
        // std::sort(std::execution::par, Z_pre.begin(), Z_pre.end(), comp);
        // std::sort(std::execution::par, W_pre.begin(), W_pre.end(), comp);

        // std::sort(std::execution::par, X_omega_new.begin(), X_omega_new.end(), comp);
        // std::sort(std::execution::par, Y_omega_new.begin(), Y_omega_new.end(), comp);
        // std::sort(std::execution::par, Z_omega_new.begin(), Z_omega_new.end(), comp);
        // std::sort(std::execution::par, W_omega_new.begin(), W_omega_new.end(), comp);
        // std::sort(std::execution::par, X_omega_pre.begin(), X_omega_pre.end(), comp);
        // std::sort(std::execution::par, Y_omega_pre.begin(), Y_omega_pre.end(), comp);
        // std::sort(std::execution::par, Z_omega_pre.begin(), Z_omega_pre.end(), comp);
        // std::sort(std::execution::par, W_omega_pre.begin(), W_omega_pre.end(), comp);
        
        std::vector<ZRoot2> XY_new(X_new_size*Y_new_size + X_new_size*Y_pre_size + X_pre_size*Y_new_size + 
                                   X_omega_new_size*Y_omega_new_size + X_omega_new_size*Y_omega_pre_size + X_omega_pre_size*Y_omega_new_size);
        std::vector<ZRoot2> ZW_new(Z_new_size*W_new_size + Z_new_size*W_pre_size + Z_pre_size*W_new_size +
                                   Z_omega_new_size*W_omega_new_size + Z_omega_new_size*W_omega_pre_size + Z_omega_pre_size*W_omega_new_size);


        // XY_newとZW_newを計算
#pragma omp parallel sections
        {
#pragma omp section
        {
            const long long shift = 0;
#pragma omp parallel for
            for(long long i = 0; i < X_new_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < Y_new_size; j++){
                    XY_new[i*Y_new_size + j + shift] = X_new[i]*X_new[i] + Y_new[j]*Y_new[j];
                }
            }
        }
#pragma omp section
        {
            const long long shift = X_new_size*Y_new_size;
#pragma omp parallel for
            for(long long i = 0; i < X_new_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < Y_pre_size; j++){
                    XY_new[i*Y_pre_size + j + shift] = X_new[i]*X_new[i] + Y_pre[j]*Y_pre[j];
                }
            }
        }
#pragma omp section
        {
            const long long shift = X_new_size*Y_new_size + X_new_size*Y_pre_size;
#pragma omp parallel for
            for(long long i = 0; i < X_pre_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < Y_new_size; j++){
                    XY_new[i*Y_new_size + j + shift] = X_pre[i]*X_pre[i] + Y_new[j]*Y_new[j];
                }
            }
        }
#pragma omp section
        {
            const long long shift = X_new_size*Y_new_size + X_new_size*Y_pre_size + X_pre_size*Y_new_size;
#pragma omp parallel for
            for(long long i = 0; i < X_omega_new_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < Y_omega_new_size; j++){
                    XY_new[i*Y_omega_new_size + j + shift] = X_omega_new[i]*X_omega_new[i] + Y_omega_new[j]*Y_omega_new[j];
                }
            }
        }
#pragma omp section
        {
            const long long shift = X_new_size*Y_new_size + X_new_size*Y_pre_size + X_pre_size*Y_new_size + 
                                    X_omega_new_size*Y_omega_new_size;
#pragma omp parallel for
            for(long long i = 0; i < X_omega_new_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < Y_omega_pre_size; j++){
                    XY_new[i*Y_omega_pre_size + j + shift] = X_omega_new[i]*X_omega_new[i] + Y_omega_pre[j]*Y_omega_pre[j];
                }
            }
        }
#pragma omp section
        {
            const long long shift = X_new_size*Y_new_size + X_new_size*Y_pre_size + X_pre_size*Y_new_size + 
                                    X_omega_new_size*Y_omega_new_size + X_omega_new_size*Y_omega_pre_size;
#pragma omp parallel for
            for(long long i = 0; i < X_omega_pre_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < Y_omega_new_size; j++){
                    XY_new[i*Y_omega_new_size + j + shift] = X_omega_pre[i]*X_omega_pre[i] + Y_omega_new[j]*Y_omega_new[j];
                }
            }
        }


#pragma omp section
        {
            const long long shift = 0;
#pragma omp parallel for
            for(long long i = 0; i < Z_new_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < W_new_size; j++){
                    ZW_new[i*W_new_size + j + shift] = Z_new[i]*Z_new[i] + W_new[j]*W_new[j];
                }
            }
        }
#pragma omp section
        {
            const long long shift = Z_new_size*W_new_size;
#pragma omp parallel for
            for(long long i = 0; i < Z_new_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < W_pre_size; j++){
                    ZW_new[i*W_pre_size + j + shift] = Z_new[i]*Z_new[i] + W_pre[j]*W_pre[j];
                }
            }
        }
#pragma omp section
        {
            const long long shift = Z_new_size*W_new_size + Z_new_size*W_pre_size;
#pragma omp parallel for
            for(long long i = 0; i < Z_pre_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < W_new_size; j++){
                    ZW_new[i*W_new_size + j + shift] = Z_pre[i]*Z_pre[i] + W_new[j]*W_new[j];
                }
            }
        }
#pragma omp section
        {
            const long long shift = Z_new_size*W_new_size + Z_new_size*W_pre_size + Z_pre_size*W_new_size;
#pragma omp parallel for
            for(long long i = 0; i < Z_omega_new_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < W_omega_new_size; j++){
                    ZW_new[i*W_omega_new_size + j + shift] = Z_omega_new[i]*Z_omega_new[i] + W_omega_new[j]*W_omega_new[j];
                }
            }
        }
#pragma omp section
        {
            const long long shift = Z_new_size*W_new_size + Z_new_size*W_pre_size + Z_pre_size*W_new_size + 
                                    Z_omega_new_size*W_omega_new_size;
#pragma omp parallel for
            for(long long i = 0; i < Z_omega_new_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < W_omega_pre_size; j++){
                    ZW_new[i*W_omega_pre_size + j + shift] = Z_omega_new[i]*Z_omega_new[i] + W_omega_pre[j]*W_omega_pre[j];
                }
            }
        }
#pragma omp section
        {
            const long long shift = Z_new_size*W_new_size + Z_new_size*W_pre_size + Z_pre_size*W_new_size + 
                                    Z_omega_new_size*W_omega_new_size + Z_omega_new_size*W_omega_pre_size;
#pragma omp parallel for
            for(long long i = 0; i < Z_omega_pre_size; i++){
#pragma omp parallel for
                for(long long j = 0; j < W_omega_new_size; j++){
                    ZW_new[i*W_omega_new_size + j + shift] = Z_omega_pre[i]*Z_omega_pre[i] + W_omega_new[j]*W_omega_new[j];
                }
            }
        }
        }


        start = std::chrono::system_clock::now();

        // XY_newとZW_newをソート
        // cout << XY_new.size() << endl;
        // cout << ZW_new.size() << endl;
#pragma omp parallel sections
        {
#pragma omp section
            std::sort(std::execution::par, XY_new.begin(), XY_new.end());
#pragma omp section
            std::sort(std::execution::par, ZW_new.begin(), ZW_new.end());
        }
        // std::sort(std::execution::par, XY_new.begin(), XY_new.end());
        // std::sort(std::execution::par, ZW_new.begin(), ZW_new.end());

        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        // cout << "ソート " << time << endl;
    
        start = std::chrono::system_clock::now();

        std::vector<ZRoot2> XY_pre(XY.size());
        std::vector<ZRoot2> ZW_pre(ZW.size());
        // cout << "XY_preサイズ " << XY_pre.size() << endl;
        // cout << "ZW_preサイズ " << ZW_pre.size() << endl;
        // XY_pre = std::move(XY);
        // ZW_pre = std::move(ZW);
        // std::copy(XY.begin(), XY.end(), XY_pre.begin());
        // std::copy(ZW.begin(), ZW.end(), ZW_pre.begin());
        XY.swap(XY_pre);
        ZW.swap(ZW_pre);

        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        // cout << "ムーブ " << time << endl;

        start = std::chrono::system_clock::now();      
        
        // 解の探索
        std::vector<ZRoot2> X(X_total.begin(), X_total.end());
        X.insert(X.end(), X_omega_total.begin(), X_omega_total.end());
        std::vector<ZRoot2> Y(Y_total.begin(), Y_total.end());
        Y.insert(Y.end(), Y_omega_total.begin(), Y_omega_total.end());
        std::vector<ZRoot2> Z(Z_total.begin(), Z_total.end());
        Z.insert(Z.end(), Z_omega_total.begin(), Z_omega_total.end());
        std::vector<ZRoot2> W(W_total.begin(), W_total.end());
        W.insert(W.end(), W_omega_total.begin(), W_omega_total.end());

        std::vector<ZRoot2> X_squared(X.size());
        for(int i = 0; i < X.size(); i++) X_squared[i] = X[i] * X[i];
        std::vector<ZRoot2> Y_squared(Y.size());
        for(int i = 0; i < Y.size(); i++) Y_squared[i] = Y[i] * Y[i];
        std::vector<ZRoot2> Z_squared(Z.size());
        for(int i = 0; i < Z.size(); i++) Z_squared[i] = Z[i] * Z[i];
        std::vector<ZRoot2> W_squared(W.size());
        for(int i = 0; i < W.size(); i++) W_squared[i] = W[i] * W[i];

        std::sort(std::execution::par, X_squared.begin(), X_squared.end());
        std::sort(std::execution::par, Y_squared.begin(), Y_squared.end());
        std::sort(std::execution::par, Z_squared.begin(), Z_squared.end());
        std::sort(std::execution::par, W_squared.begin(), W_squared.end());

        ZRoot2 key = {(ITYPE)1 << (k+1), 0};
        std::vector<std::array<ZRoot2, 4>> XYZW1, XYZW2, XYZW3;
#pragma omp parallel sections
        {
#pragma omp section
            XYZW1 = subroutine(X, Y, Z, W, X_squared, Y_squared, Z_squared, W_squared, XY_new, ZW_new, key);
#pragma omp section
            XYZW2 = subroutine(X, Y, Z, W, X_squared, Y_squared, Z_squared, W_squared, XY_pre, ZW_new, key);
#pragma omp section
            XYZW3 = subroutine(X, Y, Z, W, X_squared, Y_squared, Z_squared, W_squared, XY_new, ZW_pre, key);
        }

        // √2(x + iy) = √2(aω^3 + bω^2 + cω + d)の(a,b,c,d)を求める。
        auto convert = [](ZRoot2 re, ZRoot2 im) -> ZOmega{ return {(-re.a+im.a) / 2, im.b, (re.a+im.a) / 2, re.b}; };
        for(auto [x, y, z, w] : XYZW1) candidates.push_back({convert(x,y), convert(z,w), l, k});
        for(auto [x, y, z, w] : XYZW2) candidates.push_back({convert(x,y), convert(z,w), l, k});
        for(auto [x, y, z, w] : XYZW3) candidates.push_back({convert(x,y), convert(z,w), l, k});
        
        std::vector<U2_ZOmega> solutions;
        std::vector<U2_ZOmega> next_candidates;
        for(auto V_ZOmega : candidates){
            quaternion V = to_quaternion(V_ZOmega);
            if(distance(U, V) <= eps){
                bool flag = true;
                for(auto &sol : solutions) if(is_equal(V_ZOmega, sol)) flag = false;
                if(flag){
                    solutions.push_back(V_ZOmega);
                    continue;
                }
            }
            next_candidates.push_back(V_ZOmega);
        }
        candidates.swap(next_candidates);
        
        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        // cout << "サーチ " << time << endl;

        start = std::chrono::system_clock::now();

        XY.resize(XY_pre.size() + XY_new.size());
        ZW.resize(ZW_pre.size() + ZW_new.size());
        std::merge(std::execution::par, XY_pre.begin(), XY_pre.end(), XY_new.begin(), XY_new.end(), XY.begin());
        std::merge(std::execution::par, ZW_pre.begin(), ZW_pre.end(), ZW_new.begin(), ZW_new.end(), ZW.begin());

        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        // cout << "マージ " << time << endl;

        return solutions;
    }




    std::vector<U2_ZOmega>
    enum_u_t_use_hash
    (quaternion U, FTYPE eps, FTYPE eps_pre, int k, int l, tbb::concurrent_hash_map<ZRoot2, std::pair<ZRoot2, ZRoot2>, ZRoot2_hash>& hash_XY, tbb::concurrent_hash_map<ZRoot2, std::pair<ZRoot2, ZRoot2>, ZRoot2_hash>& hash_ZW)
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

        FTYPE eps_prime = sqrt(1 - sqrt(1 - eps * eps)) * sqrt2;
        if(isnan(eps_prime)) eps_prime = eps;

        FTYPE eps_pre_prime = sqrt(1 - sqrt(1 - eps_pre * eps_pre)) * sqrt2;
        if(isnan(eps_pre_prime)) eps_pre_prime = eps_pre;

        std::vector<ZRoot2> X_total = one_dim_grid_problem((a - eps_prime)*sqrt2k, (a + eps_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> Y_total = one_dim_grid_problem((b - eps_prime)*sqrt2k, (b + eps_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> Z_total = one_dim_grid_problem((c - eps_prime)*sqrt2k, (c + eps_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> W_total = one_dim_grid_problem((d - eps_prime)*sqrt2k, (d + eps_prime)*sqrt2k, y0, y1);
        std::sort(X_total.begin(), X_total.end());
        std::sort(Y_total.begin(), Y_total.end());
        std::sort(Z_total.begin(), Z_total.end());
        std::sort(W_total.begin(), W_total.end());

        std::vector<ZRoot2> X_pre = one_dim_grid_problem((a - eps_pre_prime)*sqrt2k, (a + eps_pre_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> Y_pre = one_dim_grid_problem((b - eps_pre_prime)*sqrt2k, (b + eps_pre_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> Z_pre = one_dim_grid_problem((c - eps_pre_prime)*sqrt2k, (c + eps_pre_prime)*sqrt2k, y0, y1);
        std::vector<ZRoot2> W_pre = one_dim_grid_problem((d - eps_pre_prime)*sqrt2k, (d + eps_pre_prime)*sqrt2k, y0, y1);
        std::sort(X_pre.begin(), X_pre.end());
        std::sort(Y_pre.begin(), Y_pre.end());
        std::sort(Z_pre.begin(), Z_pre.end());
        std::sort(W_pre.begin(), W_pre.end());

        std::vector<ZRoot2> X_new, Y_new, Z_new, W_new;
        std::set_difference(X_total.begin(), X_total.end(), X_pre.begin(), X_pre.end(), std::back_inserter(X_new));
        std::set_difference(Y_total.begin(), Y_total.end(), Y_pre.begin(), Y_pre.end(), std::back_inserter(Y_new));
        std::set_difference(Z_total.begin(), Z_total.end(), Z_pre.begin(), Z_pre.end(), std::back_inserter(Z_new));
        std::set_difference(W_total.begin(), W_total.end(), W_pre.begin(), W_pre.end(), std::back_inserter(W_new));

        y0 += inv_sqrt2;
        y1 += inv_sqrt2;
        std::vector<ZRoot2> X_omega_total = one_dim_grid_problem((a - eps_prime)*sqrt2k - inv_sqrt2, (a + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> Y_omega_total = one_dim_grid_problem((b - eps_prime)*sqrt2k - inv_sqrt2, (b + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> Z_omega_total = one_dim_grid_problem((c - eps_prime)*sqrt2k - inv_sqrt2, (c + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> W_omega_total = one_dim_grid_problem((d - eps_prime)*sqrt2k - inv_sqrt2, (d + eps_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::sort(X_omega_total.begin(), X_omega_total.end());
        std::sort(Y_omega_total.begin(), Y_omega_total.end());
        std::sort(Z_omega_total.begin(), Z_omega_total.end());
        std::sort(W_omega_total.begin(), W_omega_total.end());

        std::vector<ZRoot2> X_omega_pre = one_dim_grid_problem((a - eps_pre_prime)*sqrt2k - inv_sqrt2, (a + eps_pre_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> Y_omega_pre = one_dim_grid_problem((b - eps_pre_prime)*sqrt2k - inv_sqrt2, (b + eps_pre_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> Z_omega_pre = one_dim_grid_problem((c - eps_pre_prime)*sqrt2k - inv_sqrt2, (c + eps_pre_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::vector<ZRoot2> W_omega_pre = one_dim_grid_problem((d - eps_pre_prime)*sqrt2k - inv_sqrt2, (d + eps_pre_prime)*sqrt2k - inv_sqrt2, y0, y1);
        std::sort(X_omega_pre.begin(), X_omega_pre.end());
        std::sort(Y_omega_pre.begin(), Y_omega_pre.end());
        std::sort(Z_omega_pre.begin(), Z_omega_pre.end());
        std::sort(W_omega_pre.begin(), W_omega_pre.end());

        std::vector<ZRoot2> X_omega_new, Y_omega_new, Z_omega_new, W_omega_new;
        std::set_difference(X_omega_total.begin(), X_omega_total.end(), X_omega_pre.begin(), X_omega_pre.end(), std::back_inserter(X_omega_new));
        std::set_difference(Y_omega_total.begin(), Y_omega_total.end(), Y_omega_pre.begin(), Y_omega_pre.end(), std::back_inserter(Y_omega_new));
        std::set_difference(Z_omega_total.begin(), Z_omega_total.end(), Z_omega_pre.begin(), Z_omega_pre.end(), std::back_inserter(Z_omega_new));
        std::set_difference(W_omega_total.begin(), W_omega_total.end(), W_omega_pre.begin(), W_omega_pre.end(), std::back_inserter(W_omega_new));


        const long unsigned int X_new_size = X_new.size();
        const long unsigned int Y_new_size = Y_new.size();
        const long unsigned int Z_new_size = Z_new.size();
        const long unsigned int W_new_size = W_new.size();

        const long unsigned int X_pre_size = X_pre.size();
        const long unsigned int Y_pre_size = Y_pre.size();
        const long unsigned int Z_pre_size = Z_pre.size();
        const long unsigned int W_pre_size = W_pre.size();

        const long unsigned int X_omega_new_size = X_omega_new.size();
        const long unsigned int Y_omega_new_size = Y_omega_new.size();
        const long unsigned int Z_omega_new_size = Z_omega_new.size();
        const long unsigned int W_omega_new_size = W_omega_new.size();

        const long unsigned int X_omega_pre_size = X_omega_pre.size();
        const long unsigned int Y_omega_pre_size = Y_omega_pre.size();
        const long unsigned int Z_omega_pre_size = Z_omega_pre.size();
        const long unsigned int W_omega_pre_size = W_omega_pre.size();

        std::cout << X_new_size + X_pre_size << " " << Y_new_size + Y_pre_size << " " << Z_new_size + Z_pre_size << " " << W_new_size + W_pre_size << std::endl; 
        std::cout << X_new_size << " " << Y_new_size << " " << Z_new_size << " " << W_new_size << std::endl; 
        // std::cout << X_omega_new_size + X_omega_pre_size << " " << Y_omega_new_size + Y_omega_pre_size << " " << Z_omega_new_size + Z_omega_pre_size << " " << W_omega_new_size + W_omega_pre_size << std::endl; 

        // ((√2x)^2 + (√2x)^2 + (√2x)^2 + (√2x)^2) = 2^{k+1}を満たす(x,y,z,w)を求める
        ZRoot2 sqrt2_ZRoot2 = {0,1};
        for(auto &x : X_new) x *= sqrt2_ZRoot2;
        for(auto &y : Y_new) y *= sqrt2_ZRoot2;
        for(auto &z : Z_new) z *= sqrt2_ZRoot2;
        for(auto &w : W_new) w *= sqrt2_ZRoot2;

        for(auto &x : X_pre) x *= sqrt2_ZRoot2;
        for(auto &y : Y_pre) y *= sqrt2_ZRoot2;
        for(auto &z : Z_pre) z *= sqrt2_ZRoot2;
        for(auto &w : W_pre) w *= sqrt2_ZRoot2;

        for(auto &x : X_omega_new) x = x * sqrt2_ZRoot2 + 1;
        for(auto &y : Y_omega_new) y = y * sqrt2_ZRoot2 + 1;
        for(auto &z : Z_omega_new) z = z * sqrt2_ZRoot2 + 1;
        for(auto &w : W_omega_new) w = w * sqrt2_ZRoot2 + 1;

        for(auto &x : X_omega_pre) x = x * sqrt2_ZRoot2 + 1;
        for(auto &y : Y_omega_pre) y = y * sqrt2_ZRoot2 + 1;
        for(auto &z : Z_omega_pre) z = z * sqrt2_ZRoot2 + 1;
        for(auto &w : W_omega_pre) w = w * sqrt2_ZRoot2 + 1;

        std::vector<std::array<ZRoot2, 4>> xyzw;   // ((√2x)^2 + (√2x)^2 + (√2x)^2 + (√2x)^2) = 2^{k+1}を満たす(x,y,z,w)
        ZRoot2 R = {((ITYPE)1 << (k+1)), 0};   // R := 2^{k+1}

        /*
        (X_new, Y_new, Z_pre, W_pre), (X_new, Y_pre, Z_pre, W_pre), (X_pre, Y_pre, Z_pre, W_pre)の組み合わせを計算
        また同時にhash_XYに(X_new, Y_new), (X_new, Y_pre), (X_pre, Y_new)の組み合わせを挿入
        */

        start = std::chrono::system_clock::now();
#pragma omp parallel sections
        {
#pragma omp section
        {
            // X_newとY_new
#pragma omp parallel for
            for(auto &x : X_new){
#pragma omp parallel for
                for(auto &y : Y_new){
                    ZRoot2 x2y2 = x*x + y*y;
                    ZRoot2 key = R - x2y2;
                    // auto range = hash_ZW.equal_range(key);
                    // for(auto it = range.first; it != range.second; it++){
                    //     auto [z,w] = it->second;
                    //     if(x2y2 + z*z + w*w == R) xyzw.push_back({x,y,z,w}); 
                    // }

                    hash_XY.insert({x2y2, {x, y}});
                }
            }
        }
#pragma omp section
        {
            // X_newとY_pre
#pragma omp parallel for
            for(auto &x : X_new){
#pragma omp parallel for
                for(auto &y : Y_pre){
                    ZRoot2 x2y2 = x*x + y*y;
                    ZRoot2 key = R - x2y2;
                    auto range = hash_ZW.equal_range(key);
                    for(auto it = range.first; it != range.second; it++){
                        auto [z,w] = it->second;
                        if(x2y2 + z*z + w*w == R) xyzw.push_back({x,y,z,w}); 
                    }

                    hash_XY.insert({x2y2, {x, y}});
                }
            }
        }

#pragma omp section
        {
            // X_preとY_new
#pragma omp parallel for
            for(auto &x : X_pre){
#pragma omp parallel for
                for(auto &y : Y_new){
                    ZRoot2 x2y2 = x*x + y*y;
                    ZRoot2 key = R - x2y2;
                    auto range = hash_ZW.equal_range(key);
                    for(auto it = range.first; it != range.second; it++){
                        auto [z,w] = it->second;
                        if(x2y2 + z*z + w*w == R) xyzw.push_back({x,y,z,w}); 
                    }

                    hash_XY.insert({x2y2, {x, y}});
                }
            }
        }
#pragma omp section
        {
            // X_omega_newとY_omega_new
#pragma omp parallel for
            for(auto &x : X_omega_new){
#pragma omp parallel for
                for(auto &y : Y_omega_new){
                    ZRoot2 x2y2 = x*x + y*y;
                    ZRoot2 key = R - x2y2;
                    // auto range = hash_ZW.equal_range(key);
                    // for(auto it = range.first; it != range.second; it++){
                    //     auto [z,w] = it->second;
                    //     if(x2y2 + z*z + w*w == R) xyzw.push_back({x,y,z,w}); 
                    // }

                    hash_XY.insert({x2y2, {x, y}});
                }
            }
        }
#pragma omp section
        {
            // X_omega_newとY_omega_pre
#pragma omp parallel for
            for(auto &x : X_omega_new){
#pragma omp parallel for
                for(auto &y : Y_omega_pre){
                    ZRoot2 x2y2 = x*x + y*y;
                    ZRoot2 key = R - x2y2;
                    auto range = hash_ZW.equal_range(key);
                    for(auto it = range.first; it != range.second; it++){
                        auto [z,w] = it->second;
                        if(x2y2 + z*z + w*w == R) xyzw.push_back({x,y,z,w}); 
                    }

                    hash_XY.insert({x2y2, {x, y}});
                }
            }
        }
#pragma omp section
        {
            // X_omega_preとY_omega_new
#pragma omp parallel for
            for(auto &x : X_omega_pre){
#pragma omp parallel for
                for(auto &y : Y_omega_new){
                    ZRoot2 x2y2 = x*x + y*y;
                    ZRoot2 key = R - x2y2;
                    auto range = hash_ZW.equal_range(key);
                    for(auto it = range.first; it != range.second; it++){
                        auto [z,w] = it->second;
                        if(x2y2 + z*z + w*w == R) xyzw.push_back({x,y,z,w}); 
                    }

                    hash_XY.insert({x2y2, {x, y}});
                }
            }
        }
        }
        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        std::cout << "XYの時間 " << time << "[ms]" << std::endl; 


        start = std::chrono::system_clock::now();
        /*
        (X_total, Y_total, Z_new, W_new), (X_total, Y_total, Z_new, W_pre), (X_total, Y_total, Z_pre, W_pre)の組み合わせを計算
        また同時にhash_ZWに(Z_new, W_new), (Z_new, W_pre), (Z_pre, W_new)の組み合わせを挿入
        */
#pragma omp parallel sections
        {
#pragma omp section
        {
            // Z_newとW_new
            for(auto &z : Z_new){
                for(auto &w : W_new){
                    ZRoot2 z2w2 = z*z + w*w;
                    ZRoot2 key = R - z2w2;
                    auto range = hash_XY.equal_range(key);
                    for(auto it = range.first; it != range.second; it++){
                        auto [x,y] = it->second;
                        if(x*x + y*y + z2w2 == R) xyzw.push_back({x,y,z,w});
                    }

                    hash_ZW.insert({z2w2, {z, w}});
                }
            }
        }
#pragma omp section
        {
            // Z_newとW_pre
            for(auto &z : Z_new){
                for(auto &w : W_pre){
                    ZRoot2 z2w2 = z*z + w*w;
                    ZRoot2 key = R - z2w2;
                    auto range = hash_XY.equal_range(key);
                    for(auto it = range.first; it != range.second; it++){
                        auto [x,y] = it->second;
                        if(x*x + y*y + z2w2 == R) xyzw.push_back({x,y,z,w});
                    }

                    hash_ZW.insert({z2w2, {z, w}});
                }
            }
        }
#pragma omp section
        {
            // Z_preとW_new
            for(auto &z : Z_pre){
                for(auto &w : W_new){
                    ZRoot2 z2w2 = z*z + w*w;
                    ZRoot2 key = R - z2w2;
                    auto range = hash_XY.equal_range(key);
                    for(auto it = range.first; it != range.second; it++){
                        auto [x,y] = it->second;
                        if(x*x + y*y + z2w2 == R) xyzw.push_back({x,y,z,w});
                    }

                    hash_ZW.insert({z2w2, {z, w}});
                }
            }
        }
#pragma omp section
        {
            // Z_omega_newとW_omega_new
            for(auto &z : Z_omega_new){
                for(auto &w : W_omega_new){
                    ZRoot2 z2w2 = z*z + w*w;
                    ZRoot2 key = R - z2w2;
                    auto range = hash_XY.equal_range(key);
                    for(auto it = range.first; it != range.second; it++){
                        auto [x,y] = it->second;
                        if(x*x + y*y + z2w2 == R) xyzw.push_back({x,y,z,w});
                    }

                    hash_ZW.insert({z2w2, {z, w}});
                }
            }
        }
#pragma omp section
        {
            // Z_omega_newとW_omega_pre
            for(auto &z : Z_omega_new){
                for(auto &w : W_omega_pre){
                    ZRoot2 z2w2 = z*z + w*w;
                    ZRoot2 key = R - z2w2;
                    auto range = hash_XY.equal_range(key);
                    for(auto it = range.first; it != range.second; it++){
                        auto [x,y] = it->second;
                        if(x*x + y*y + z2w2 == R) xyzw.push_back({x,y,z,w});
                    }

                    hash_ZW.insert({z2w2, {z, w}});
                }
            }
        }
#pragma omp section
        {
            // Z_omega_preとW_omega_new
            for(auto &z : Z_omega_pre){
                for(auto &w : W_omega_new){
                    ZRoot2 z2w2 = z*z + w*w;
                    ZRoot2 key = R - z2w2;
                    auto range = hash_XY.equal_range(key);
                    for(auto it = range.first; it != range.second; it++){
                        auto [x,y] = it->second;
                        if(x*x + y*y + z2w2 == R) xyzw.push_back({x,y,z,w});
                    }

                    hash_ZW.insert({z2w2, {z, w}});
                }
            }
        }
        }
        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        std::cout << "ZWの時間 " << time << "[ms]" << std::endl;


        // u = √2(x+iy)をZ[ω]の係数に変換する関数.
        auto convert = [](ZRoot2 re, ZRoot2 im) -> ZOmega{ return {(-re.a+im.a) / 2, im.b, (re.a+im.a) / 2, re.b}; };


        std::vector<U2_ZOmega> candidates;
        for(auto [x,y,z,w] : xyzw) candidates.push_back({convert(x, y), convert(z, w), l, k});

        std::cout << "候補数 " << candidates.size() << std::endl;

        std::vector<U2_ZOmega> ret;
        for(auto &V_ZOmega : candidates){
            quaternion V = to_quaternion(V_ZOmega);
            if(distance(U, V) <= eps) ret.push_back(V_ZOmega);
        }

        return ret;
    }

}