#include <bits/stdc++.h>

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
    b << std::sqrt(1 - (eps*eps / 4)), std::sqrt(1 - (eps*eps / 4));
    
    auto[Ar, br] = row_reduction<T, 2, 4>(A, b);

    Eigen::Matrix<T, 4, 1> inter;   // Ax=bが張る平面の切片
    inter = Eigen::Matrix<T, 4, 1>::Zero(4);
    for(int i = 0; i < 2; i++) inter(i) = br(i);   // intercept=[br1, br2, 0, 0]
    Eigen::Matrix<T, 4, 2> V = Gram_Schmidt<T, 4, 2>(A.transpose());   // Ax=0が張る平面の直交補空間の正規直交基底

    T v0_dot = inter.dot(V.col(0));
    T v1_dot = inter.dot(V.col(1));
    T dist = std::sqrt(v0_dot*v0_dot + v1_dot*v1_dot);   // Ax=bが張る平面と原点との距離

    if(dist <= (T)1.0) return true;
    else return false;
}


/*
B(u1,ε)とB(u2,ε)が交わるある１点を計算する。
出力 : 交点を１つを出力する
*/
template<typename T>
bool call_inter_2ball(quaternion<T> u1, quaternion<T> u2, T eps)
{
    Eigen::Matrix<T, 2, 4> A;
    Eigen::Matrix<T, 2, 1> b;
    auto [a11, a12, a13, a14] = u1.get_abcd();
    auto [a21, a22, a23, a24] = u2.get_abcd();
    A << a11, a12, a13, a14,
         a21, a22, a23, a24;
    b << std::sqrt(1 - (eps*eps / 4)), std::sqrt(1 - (eps*eps / 4));
    
    auto[Ar, br] = row_reduction<T, 2, 4>(A, b);

    Eigen::Matrix<T, 4, 1> v;
    Eigen::Matrix<T, 4, 1> inter;   // Ax=bが張る平面の切片
    inter = Eigen::Matrix<T, 4, 1>::Zero(4);

    for(int i = 0; i < 2; i++) v(i) = -Ar(i, 3);
    for(int i = 0; i < 2; i++) inter(i) = br(i);   // intercept=[br1, br2, 0, 0]

    T v0_dot = inter.dot(V.col(0));
    T v1_dot = inter.dot(V.col(1));
    T dist = std::sqrt(v0_dot*v0_dot + v1_dot*v1_dot);   // Ax=bが張る平面と原点との距離

    if(dist <= (T)1.0) return true;
    else return false;
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
    b << std::sqrt(1 - (eps*eps / 4)), std::sqrt(1 - (eps*eps / 4)), std::sqrt(1 - (eps*eps / 4));
    
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

    T t1 = (-bb + std::sqrt(D)) / (2*a);
    T t2 = (-bb - std::sqrt(D)) / (2*a);

    Eigen::Matrix<T, 4, 1> ret1_eigen, ret2_eigen;
    ret1_eigen = t1*v + inter;
    ret2_eigen = t2*v + inter;

    quaternion<T> ret1 = convert_quaternion(ret1_eigen);
    quaternion<T> ret2 = convert_quaternion(ret2_eigen);

    return {ret1, ret2};
}


template<typename T>
bool check_eps_net(std::vector<quaternion<T>> X, std::vector<T> U_target, T eps){  
    int N = X.szie();
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i+1; j++){
            if(!check_inter_2ball(X[i], X[j], eps)) continue;
            bool flag_cross = false;
            for(int k = j+1; k < N; k++){
                auto [inter1, inter2] = cal_inter_3ball(X(i), X(j), X(k));
                if(!inter1.is_unitrary() || !inter2.is_unitary()) continue;
                
                if(distance(inter1, U_target) < 2*eps){   // 2ε-ballに含まれる
                    flag_cross = true;
                    bool flag_in = false;
                    for(int l = 0; l < N; l++){
                        if(l == i || l == j || l == k) continue;
                        if(distance(inter1, X[l])) flag_in = true; 
                    }
                    if(flag_in == false) return false;
                }

                if(distance(inter2, U_target) < 2*eps){   // 2ε-ballに含まれる
                    flag_cross = true;
                    bool flag_in = false;
                    for(int l = 0; l < N; l++){
                        if(l == i || l == j || l == k) continue;
                        if(distance(inter2, X[l])) flag_in = true; 
                    }
                    if(flag_in == false) return false;
                }
            }   // kのループ
            if(!flag_cross){
                bool flag_in = false;
                quaternion<T> 
            }
        }
    }
}