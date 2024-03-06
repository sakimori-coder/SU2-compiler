#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "quaternion.hpp"
#include "linalg.hpp"
#include <Eigen/Core>

using std::cout;
using std::endl;


/*
B(u1,ε)とB(u2,ε)が交わるかを判定する。
出力 : 交わるならtrue, 交わらないならfalse
*/
template<typename T>
bool check_inter_2ball(quaternion<T> u1, quaternion<T> u2, T eps)
{
    Eigen::Matrix<T, 2, 4> A;
    Eigen::Matrix<T, 2, 1> b;
    auto [a11, a12, a13, a14] = u1.get_abcd();
    auto [a21, a22, a23, a24] = u2.get_abcd();
    A << a11, a12, a13, a14,
         a21, a22, a23, a24;
    T tmp = (T)1.0 - eps*eps;
    b << sqrt(tmp), sqrt(tmp);

    auto [V, inter] = solve_linear_system<T, 2, 4>(A, b);

    Eigen::Matrix<T, 4, 2> OC = Gram_Schmidt<T, 4, 2>(A.transpose());   // Ax=0が張る平面の直交補空間の正規直交基底

    T dot0 = inter.dot(OC.col(0));
    T dot1 = inter.dot(OC.col(1));
    T dist = sqrt(dot0*dot0 + dot1*dot1);   // Ax=bが張る平面と原点との距離
    
    if(dist <= (T)1.0) return true;
    else return false;
}


/*
B(u1,ε)とB(u2,ε)が交わるある１点を計算する。
出力 : 交点を１つを出力する
*/
template<typename T>
quaternion<T> cal_inter_2ball(quaternion<T> u1, quaternion<T> u2, T eps)
{
    Eigen::Matrix<T, 2, 4> A;
    Eigen::Matrix<T, 2, 1> b;
    auto [a11, a12, a13, a14] = u1.get_abcd();
    auto [a21, a22, a23, a24] = u2.get_abcd();
    A << a11, a12, a13, a14,
         a21, a22, a23, a24;
    T tmp = (T)1.0 - eps*eps;
    b << sqrt(tmp), sqrt(tmp);

    auto [V, inter] = solve_linear_system<T, 2, 4>(A, b);

    Eigen::Matrix<T, 4, 1> v1, v2;
    v1 = V.col(0);
    v2 = V.col(1);

    Eigen::Matrix<T, 2, 1> mu;
    mu << v1.dot(inter), v2.dot(inter);
    T s = sqrt(1 - inter.norm()*inter.norm() + mu.norm()*mu.norm()) - v1.dot(inter);
    T t = -v2.dot(inter);

    // T s = -v1.dot(inter);
    // T t = sqrt(1 - inter.norm()*inter.norm() + mu.norm()*mu.norm()) - v2.dot(inter);

    Eigen::Matrix<T, 4, 1> ret_eigen;
    ret_eigen = s*v1 + t*v2 + inter;
    
    quaternion<T> ret = convert_quaternion(ret_eigen);
    
    return ret;
}


template<typename T>
quaternion<T> cal_inter_2ball(quaternion<T> u1, quaternion<T> u2, quaternion<T> u_target, T eps)
{
    Eigen::Matrix<T, 4, 1> u;
    u << u_target.get_a(), u_target.get_b(), u_target.get_c(), u_target.get_d();
    Eigen::Matrix<T, 2, 4> A;
    Eigen::Matrix<T, 2, 1> b;
    auto [a11, a12, a13, a14] = u1.get_abcd();
    auto [a21, a22, a23, a24] = u2.get_abcd();
    A << a11, a12, a13, a14,
         a21, a22, a23, a24;
    T tmp = (T)1.0 - eps*eps;
    b << sqrt(tmp), sqrt(tmp);

    auto [V, inter] = solve_linear_system<T, 2, 4>(A, b);

    Eigen::Matrix<T, 4, 1> v1, v2;
    v1 = V.col(0);
    v2 = V.col(1);

    Eigen::Matrix<T, 2, 1> mu;
    mu << v1.dot(inter), v2.dot(inter);

    T c = (T)1.0 - inter.dot(inter) + mu.dot(mu);
    inter = inter - mu(0)*v1 - mu(1)*v2;

    T u_v1 = u.dot(v1);
    T square_u_v1 = u_v1*u_v1;
    T square_u_v2 = u.dot(v2)*u.dot(v2);
    T s = sqrt(square_u_v1*c / (square_u_v1 + square_u_v2));
    T t = sqrt(c - s*s);

    std::vector<quaternion<T>> cand;
    cand.push_back(convert_quaternion<T>( s*v1 + t*v2 + inter));
    cand.push_back(convert_quaternion<T>(-s*v1 - t*v2 + inter));
    cand.push_back(convert_quaternion<T>( s*v1 - t*v2 + inter));
    cand.push_back(convert_quaternion<T>(-s*v1 + t*v2 + inter));


    std::vector<T> dist(4);
    for(int i = 0; i < 4; i++) dist[i] = distance(u_target, cand[i]);

    auto max_it = std::max_element(dist.begin(), dist.end());
    auto min_it = std::min_element(dist.begin(), dist.end());


    quaternion<T> ret = cand[0];
    T min_dist = dist[0];
    for(int i = 1; i < 4; i++){
        if(dist[i] < min_dist){
            min_dist = dist[i];
            ret = cand[i];
        }
    }

    return ret;
}



/*
B(u1,ε)とB(u2,ε)とB(u3,ε)が交わる点を計算する。
出力 : 交わるならば、交わる2点を出力する。交わらなければ,(0,0,0,0)を返す。
*/
template<typename T>
std::pair<quaternion<T>, quaternion<T>> cal_inter_3ball(quaternion<T> u1, quaternion<T> u2, quaternion<T> u3, T eps)
{
    Eigen::Matrix<T, 3, 4> A;
    Eigen::Matrix<T, 3, 1> b;
    auto [a11, a12, a13, a14] = u1.get_abcd();
    auto [a21, a22, a23, a24] = u2.get_abcd();
    auto [a31, a32, a33, a34] = u3.get_abcd();
    A << a11, a12, a13, a14, 
         a21, a22, a23, a24, 
         a31, a32, a33, a34;
    T tmp = (T)1.0 - eps*eps;
    b << sqrt(tmp), sqrt(tmp), sqrt(tmp);
    
    auto [v, inter] = solve_linear_system<T, 3, 4>(A, b);

    // 二次方程式の係数(bbは変数名が被るから)
    T a = v.norm() * v.norm();
    T bb = 2.0 * v.dot(inter);
    T c = inter.norm() * inter.norm() - 1;
    T D = bb*bb - 4*a*c;
    if(D < 0) return {{0,0,0,0}, {0,0,0,0}};

    T t1 = (-bb + sqrt(D)) / (2*a);
    T t2 = (-bb - sqrt(D)) / (2*a);

    Eigen::Matrix<T, 4, 1> ret1_eigen, ret2_eigen;
    ret1_eigen = t1*v + inter;
    ret2_eigen = t2*v + inter;

    quaternion<T> ret1 = convert_quaternion(ret1_eigen);
    quaternion<T> ret2 = convert_quaternion(ret2_eigen);

    return {ret1, ret2};
}



template<typename T>
bool check_eps_net(std::vector<quaternion<T>> X, quaternion<T> U_target, T eps){  
    int N = X.size();
    if(N == 0) return false;

    bool success = true;
#pragma omp parallel for
    for(int i = 0; i < N; i++){
        // cout << i << endl;
        bool flag_cross_i = false;   // B(x_i)がB(u_target)内で他のε-ballと交点を持つならtrue
#pragma omp parallel for
        for(int j = 0; j < N; j++){
            if(i == j) continue;
            if(!check_inter_2ball(X[i], X[j], eps)) continue;
            // この時点では分からない flag_cross_i = true;
            bool flag_cross_j = false;   // B(x_i)とB(x_j)の交じわる領域がB(u_target)内で他のε-ballと交点を持つならtrue
#pragma omp parallel for
            for(int k = 0; k < N; k++){
                if(k == i || k == j) continue;
                auto [inter1, inter2] = cal_inter_3ball(X[i], X[j], X[k], eps);

                if(!inter1.is_unitary() || !inter2.is_unitary()) continue;

                // if(distance(inter1, U_target) < 2*eps){   // 2ε-ballに含まれる
                if(distance(inter1, U_target) < eps){   // ε-ballに含まれる
                    flag_cross_i = true;
                    flag_cross_j = true;
                    bool flag_in = false;
#pragma omp parallel for
                    for(int l = 0; l < N; l++){
                        if(l == i || l == j || l == k) continue;
                        if(distance(inter1, X[l]) < eps) flag_in = true;
                    }
                    if(!flag_in) success = false;
                }

                // if(distance(inter2, U_target) < 2*eps){   // 2ε-ballに含まれる
                if(distance(inter2, U_target) < eps){   // ε-ballに含まれる
                    flag_cross_i = true;
                    flag_cross_j = true;
                    bool flag_in = false;
#pragma omp parallel for
                    for(int l = 0; l < N; l++){
                        if(l == i || l == j || l == k) continue;
                        // cout << eps - distance(inter2, X[l]) << endl; 
                        if(distance(inter2, X[l]) < eps) flag_in = true; 
                    }
                    if(!flag_in) {success = false;}
                }
            }   // kのループ

            if(!flag_cross_j){
                bool flag_in = false;
                quaternion<T> inter = cal_inter_2ball(X[i], X[j], U_target, eps);
                if(distance(inter, U_target) > eps) continue;   // ε-ballに含まれない
                flag_cross_i = true;  
#pragma omp parallel for
                for(int k = 0; k < N; k++){
                    if(k == i || k == j) continue;
                    if(distance(inter, X[k]) < eps) flag_in = true;
                }
                if(!flag_in) {success = false;} 
            }
        }   // jのループ
        if(!flag_cross_i){
            quaternion<T> inter = cal_inter_2ball(X[i], U_target, eps);
            if(!inter.is_unitary()) {success = false;}   // 一応書いてるけど、必ず交点は持つはず

            bool flag_in = false;
            for(int j = 0; j < N; j++){
                if(j == i) continue;
                if(distance(inter, X[j]) < eps) flag_in = true;
            }
            if(!flag_in) {success = false;}
            
        }
    }   // iのループ
    
    if(success) return true;
    else return false;
}




template<typename T>
bool check_eps_net_single_theread(std::vector<quaternion<T>> X, quaternion<T> U_target, T eps){  
    int N = X.size();

    for(int i = 0; i < N; i++){
        cout << i << endl;
        bool flag_cross_i = false;
        for(int j = 0; j < N; j++){
            if(i == j) continue;
            if(!check_inter_2ball(X[i], X[j], eps)) continue;
            flag_cross_i = true;
            bool flag_cross_j = false;
            for(int k = 0; k < N; k++){
                if(k == i || k == j) continue;
                auto [inter1, inter2] = cal_inter_3ball(X[i], X[j], X[k], eps);

                if(!inter1.is_unitary() || !inter2.is_unitary()) continue;
                
                if(distance(inter1, U_target) < 2*eps){   // 2ε-ballに含まれる
                    flag_cross_j = true;
                    bool flag_in = false;
                    for(int l = 0; l < N; l++){
                        if(l == i || l == j || l == k) continue;
                        if(distance(inter1, X[l]) < eps) flag_in = true;
                    }
                    if(!flag_in) return false;
                }

                if(distance(inter2, U_target) < 2*eps){   // 2ε-ballに含まれる
                    flag_cross_j = true;
                    bool flag_in = false;
                    for(int l = 0; l < N; l++){
                        if(l == i || l == j || l == k) continue;
                        if(distance(inter2, X[l]) < eps) flag_in = true; 
                    }
                    if(!flag_in) return false;
                }
            }   // kのループ

            if(!flag_cross_j){
                bool flag_in = false;
                quaternion<T> inter = cal_inter_2ball(X[i], X[j], eps);
                for(int k = 0; k < N; k++){
                    if(k == i || k == j) continue;
                    if(distance(inter, X[k]) < eps) flag_in = true;
                }
                if(!flag_in) return false; 
            }
        }   // jのループ
        if(!flag_cross_i) return false;
    }   // iのループ
    return true;
}