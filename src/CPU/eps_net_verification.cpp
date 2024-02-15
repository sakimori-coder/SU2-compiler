#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "quaternion.cpp"
#include <Eigen/Core>

using std::cout;
using std::endl;

// Ax=bの掃き出し法　計算量O(NM*min(N,M))
template<typename T, int N, int M>
std::pair<Eigen::Matrix<T, N, M>, Eigen::Matrix<T, N, 1>> 
row_reduction(Eigen::Matrix<T, N, M> A, Eigen::Matrix<T, N, 1> b)
{
    for(int k = 0; k < std::min(N, M); k++){
        T a_kk = A(k,k);
        for(int j = 0; j < M; j++) A(k,j) /= a_kk;
        b(k) /= a_kk;

        for(int i = 0; i < N; i++){
            if(i == k) continue;
            T a_ik = A(i,k);
            for(int j = 0; j < M; j++) A(i,j) -= a_ik*A(k, j);
            b(i) -= a_ik*b(k);
        }
    }

    return {A, b};
}


/*
グラム・シュミットの正規直交化  
入力 : a_i (i=1,...,M)を正規直交化したいN次元のベクトルとすると、A = (a_1,...,a_M)で与える
出力 : 正規直交化された基底を並べた行列 V = (v_1,...,v_M)を返す
*/
template<typename T, int N, int M>
Eigen::Matrix<T, N, M> Gram_Schmidt(Eigen::Matrix<T, N, M> A)
{
    Eigen::Matrix<T, N, M> V;
    for(int i = 0; i < M; i++){
        Eigen::Matrix<T, N, 1> v = A.col(i);
        for(int j = 0; j < i; j++) v -= V.col(j).dot(A.col(i)) * V.col(j);
        v = v / v.norm();
        V.col(i) = v;
    }

    return V;
}


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
    b << sqrt(1 - (eps*eps / 4)), sqrt(1 - (eps*eps / 4));
    
    auto[Ar, br] = row_reduction<T, 2, 4>(A, b);

    Eigen::Matrix<T, 4, 1> inter;   // Ax=bが張る平面の切片
    inter = Eigen::Matrix<T, 4, 1>::Zero(4);
    for(int i = 0; i < 2; i++) inter(i) = br(i);   // intercept=[br1, br2, 0, 0]
    Eigen::Matrix<T, 4, 2> V = Gram_Schmidt<T, 4, 2>(A.transpose());   // Ax=0が張る平面の直交補空間の正規直交基底

    T v0_dot = inter.dot(V.col(0));
    T v1_dot = inter.dot(V.col(1));
    T dist = sqrt(v0_dot*v0_dot + v1_dot*v1_dot);   // Ax=bが張る平面と原点との距離

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
    b << sqrt(1 - (eps*eps / 4)), sqrt(1 - (eps*eps / 4));
    
    auto[Ar, br] = row_reduction<T, 2, 4>(A, b);

    Eigen::Matrix<T, 4, 1> v1, v2;
    Eigen::Matrix<T, 4, 1> inter;   // Ax=bが張る平面の切片
    inter = Eigen::Matrix<T, 4, 1>::Zero(4);

    for(int i = 0; i < 2; i++) v1(i) = -Ar(i, 2);
    v1(2) = 1; v1(3) = 0;
    for(int i = 0; i < 2; i++) v2(i) = -Ar(i, 3);
    v2(2) = 0; v2(3) = 1;
    for(int i = 0; i < 2; i++) inter(i) = br(i);   // intercept=[br1, br2, 0, 0]

    v1 = v1 / v1.norm();
    v2 = v2 - v2.dot(v1)*v1;   // 直交化
    v2 = v2 / v2.norm();

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
    b << sqrt(1 - (eps*eps / 4)), sqrt(1 - (eps*eps / 4)), sqrt(1 - (eps*eps / 4));
    
    auto [Ar, br] = row_reduction<T, 3, 4>(A, b);
    
    Eigen::Matrix<T, 4, 1> v;       // Ax=bが張る直線の方向ベクトル
    Eigen::Matrix<T, 4, 1> inter;   // Ax=bが張る直線の切片
    inter = Eigen::Matrix<T, 4, 1>::Zero(4);
    
    for(int i = 0; i < 3; i++) v(i) = -Ar(i, 3);   // v = [-Ar03,-Ar13,-Ar23, 1]   
    v(3) = 1;
    for(int i = 0; i < 3; i++) inter(i) = br(i);   // intercept=[br1, br2, br3, 0]

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

    bool success = true;
#pragma omp parallel for
    for(int i = 0; i < N; i++){
        // cout << i << endl;
        bool flag_cross_i = false;
#pragma omp parallel for
        for(int j = 0; j < N; j++){
            if(i == j) continue;
            if(!check_inter_2ball(X[i], X[j], eps)) continue;
            flag_cross_i = true;
            bool flag_cross_j = false;
#pragma omp parallel for
            for(int k = 0; k < N; k++){
                if(k == i || k == j) continue;
                auto [inter1, inter2] = cal_inter_3ball(X[i], X[j], X[k], eps);

                if(!inter1.is_unitary() || !inter2.is_unitary()) continue;
    
                
                if(distance(inter1, U_target) < 2*eps){   // 2ε-ballに含まれる
                    flag_cross_j = true;
                    bool flag_in = false;
#pragma omp parallel for
                    for(int l = 0; l < N; l++){
                        if(l == i || l == j || l == k) continue;
                        if(distance(inter1, X[l]) < eps) flag_in = true;
                    }
                    if(!flag_in) success = false;
                }

                if(distance(inter2, U_target) < 2*eps){   // 2ε-ballに含まれる
                    flag_cross_j = true;
                    bool flag_in = false;
#pragma omp parallel for
                    for(int l = 0; l < N; l++){
                        if(l == i || l == j || l == k) continue;
                        if(distance(inter2, X[l]) < eps) flag_in = true; 
                    }
                    if(!flag_in) success = false;
                }
            }   // kのループ

            if(!flag_cross_j){
                bool flag_in = false;
                quaternion<T> inter = cal_inter_2ball(X[i], X[j], eps);
#pragma omp parallel for
                for(int k = 0; k < N; k++){
                    if(k == i || k == j) continue;
                    if(distance(inter, X[k]) < eps) flag_in = true;
                }
                if(!flag_in) success = false; 
            }
        }   // jのループ
        if(!flag_cross_i) success = false;
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