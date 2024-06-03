#include "eps_net_verification.hpp"

#include <iostream>
#include <vector>
#include <utility>
#include <Eigen/Core>

#include "type.hpp"
#include "quaternion.hpp"
#include "linalg.hpp"

using std::cout;
using std::endl;



namespace SU2_Compiler
{
    bool check_inter_2ball(quaternion U1, quaternion U2, FTYPE eps)
    {
        Eigen::Matrix<FTYPE, 2, 4> A;
        Eigen::Matrix<FTYPE, 2, 1> b;
        A << U1.a, U1.b, U1.c, U1.d,
             U2.a, U2.b, U2.c, U2.d;
        FTYPE tmp = (FTYPE)1.0 - eps*eps;
        b << sqrt(tmp), sqrt(tmp);

        auto [V, inter] = solve_linear_system<FTYPE, 2, 4>(A, b);

        Eigen::Matrix<FTYPE, 4, 2> OC = Gram_Schmidt<FTYPE, 4, 2>(A.transpose());   // Ax=0が張る平面の直交補空間の正規直交基底

        FTYPE dot0 = inter.dot(OC.col(0));
        FTYPE dot1 = inter.dot(OC.col(1));
        FTYPE dist = sqrt(dot0*dot0 + dot1*dot1);   // Ax=bが張る平面と原点との距離
        
        if(dist <= (FTYPE)1.0) return true;
        else return false;
    }


    quaternion cal_inter_2ball(quaternion U1, quaternion U2, FTYPE eps)
    {
        Eigen::Matrix<FTYPE, 2, 4> A;
        Eigen::Matrix<FTYPE, 2, 1> b;
        A << U1.a, U1.b, U1.c, U1.d,
             U2.a, U2.b, U2.c, U2.d;

        FTYPE tmp = (FTYPE)1.0 - eps*eps;
        b << sqrt(tmp), sqrt(tmp);

        auto [V, inter] = solve_linear_system<FTYPE, 2, 4>(A, b);

        Eigen::Matrix<FTYPE, 4, 1> v1, v2;
        v1 = V.col(0);
        v2 = V.col(1);

        Eigen::Matrix<FTYPE, 2, 1> mu;
        mu << v1.dot(inter), v2.dot(inter);
        FTYPE s = sqrt(1 - inter.norm()*inter.norm() + mu.norm()*mu.norm()) - v1.dot(inter);
        FTYPE t = -v2.dot(inter);

        Eigen::Matrix<FTYPE, 2, 2> AA;
        AA << v1.dot(v1), 0,
              0, v2.dot(v2);
        Eigen::Matrix<FTYPE, 2, 1> bb;
        bb = mu;
        FTYPE cc = inter.dot(inter) - 1;
        
        Eigen::Matrix<FTYPE, 2, 2> AA_inv;
        AA_inv << 1 / v1.dot(v1), 0,
                  0, 1 / v2.dot(v2);
        FTYPE cc_prime = cc - bb.dot(AA_inv * bb);


        Eigen::Matrix<FTYPE, 4, 1> ret_eigen;
        ret_eigen = s*v1 + t*v2 + inter;
        
        quaternion ret(ret_eigen(0), ret_eigen(1), ret_eigen(2), ret_eigen(3));
        
        return ret;
    }


    quaternion cal_inter_2ball(quaternion U1, quaternion U2, quaternion targetU, FTYPE eps)
    {
        Eigen::Matrix<FTYPE, 4, 1> u;
        u << targetU.a, targetU.b, targetU.c, targetU.d;
        Eigen::Matrix<FTYPE, 2, 4> A;
        Eigen::Matrix<FTYPE, 2, 1> b;
        A << U1.a, U1.b, U1.c, U1.d,
             U2.a, U2.b, U2.c, U2.d;

        FTYPE tmp = (FTYPE)1.0 - eps*eps;
        b << sqrt(tmp), sqrt(tmp);

        auto [V, inter] = solve_linear_system<FTYPE, 2, 4>(A, b);

        Eigen::Matrix<FTYPE, 4, 1> v1, v2;
        v1 = V.col(0);
        v2 = V.col(1);

        Eigen::Matrix<FTYPE, 2, 1> mu;
        mu << v1.dot(inter), v2.dot(inter);

        FTYPE c = (FTYPE)1.0 - inter.dot(inter) + mu.dot(mu);
        inter = inter - mu(0)*v1 - mu(1)*v2;

        FTYPE u_v1 = u.dot(v1);
        FTYPE square_u_v1 = u_v1*u_v1;
        FTYPE square_u_v2 = u.dot(v2)*u.dot(v2);
        FTYPE s = sqrt(square_u_v1*c / (square_u_v1 + square_u_v2));
        FTYPE t = sqrt(c - s*s);

        auto to_quaternion = [](Eigen::Matrix<FTYPE, 4, 1> x){ return quaternion(x(0), x(1), x(2), x(3)); };
        std::vector<quaternion> cand;
        cand.push_back(to_quaternion( s*v1 + t*v2 + inter));
        cand.push_back(to_quaternion(-s*v1 - t*v2 + inter));
        cand.push_back(to_quaternion( s*v1 - t*v2 + inter));
        cand.push_back(to_quaternion(-s*v1 + t*v2 + inter));

        std::vector<FTYPE> dist(4);
        for(int i = 0; i < 4; i++) dist[i] = distance(targetU, cand[i]);

        quaternion ret = cand[0];
        FTYPE min_dist = dist[0];
        for(int i = 1; i < 4; i++){
            if(dist[i] < min_dist){
                min_dist = dist[i];
                ret = cand[i];
            }
        }

        return ret;
    }


    std::pair<quaternion, quaternion> cal_inter_3ball(quaternion U1, quaternion U2, quaternion U3, FTYPE eps)
    {
        Eigen::Matrix<FTYPE, 3, 4> A;
        Eigen::Matrix<FTYPE, 3, 1> b;
        A << U1.a, U1.b, U1.c, U1.d,
             U2.a, U2.b, U2.c, U2.d,
             U3.a, U3.b, U3.c, U3.d;

        FTYPE tmp = (FTYPE)1.0 - eps*eps;
        b << sqrt(tmp), sqrt(tmp), sqrt(tmp);
        
        auto [v, inter] = solve_linear_system<FTYPE, 3, 4>(A, b);

        // 二次方程式の係数(bbは変数名が被るから)
        FTYPE a = v.norm() * v.norm();
        FTYPE bb = 2.0 * v.dot(inter);
        FTYPE c = inter.norm() * inter.norm() - 1;
        FTYPE D = bb*bb - 4*a*c;
        if(D < 0) return {{0,0,0,0}, {0,0,0,0}};

        FTYPE t1 = (-bb + sqrt(D)) / (2*a);
        FTYPE t2 = (-bb - sqrt(D)) / (2*a);

        Eigen::Matrix<FTYPE, 4, 1> ret1_eigen, ret2_eigen;
        ret1_eigen = t1*v + inter;
        ret2_eigen = t2*v + inter;

    auto to_quaternion = [](Eigen::Matrix<FTYPE, 4, 1> x){ return quaternion(x(0), x(1), x(2), x(3)); };
        quaternion ret1 = to_quaternion(ret1_eigen);
        quaternion ret2 = to_quaternion(ret2_eigen);

        return {ret1, ret2};
    }


    bool check_eps_net(std::vector<quaternion> availableU, quaternion targetU, FTYPE eps){  
        int N = availableU.size();
        if(N == 0) return false;

        bool success = true;
    #pragma omp parallel for
        for(int i = 0; i < N; i++){
            // cout << i << endl;
            bool flag_cross_i = false;   // B(x_i)がB(u_target)内で他のε-ballと交点を持つならtrue
    #pragma omp parallel for
            for(int j = 0; j < N; j++){
                if(i == j) continue;
                if(!check_inter_2ball(availableU[i], availableU[j], eps)) continue;
                // この時点では分からない flag_cross_i = true;
                bool flag_cross_j = false;   // B(x_i)とB(x_j)の交じわる領域がB(u_target)内で他のε-ballと交点を持つならtrue
    #pragma omp parallel for
                for(int k = 0; k < N; k++){
                    if(k == i || k == j) continue;
                    auto [inter1, inter2] = cal_inter_3ball(availableU[i], availableU[j], availableU[k], eps);

                    if(!inter1.is_unitary() || !inter2.is_unitary()) continue;

                    // if(distance(inter1, U_target) < 2*eps){   // 2ε-ballに含まれる
                    if(distance(inter1, targetU) < eps){   // ε-ballに含まれる
                        flag_cross_i = true;
                        flag_cross_j = true;
                        bool flag_in = false;
    #pragma omp parallel for
                        for(int l = 0; l < N; l++){
                            if(l == i || l == j || l == k) continue;
                            if(distance(inter1, availableU[l]) < eps) flag_in = true;
                        }
                        if(!flag_in){
                            success = false;
                            // std::cout << "パターン1" << std::endl;
                        }
                    }

                    // if(distance(inter2, U_target) < 2*eps){   // 2ε-ballに含まれる
                    if(distance(inter2, targetU) < eps){   // ε-ballに含まれる
                        flag_cross_i = true;
                        flag_cross_j = true;
                        bool flag_in = false;
    #pragma omp parallel for
                        for(int l = 0; l < N; l++){
                            if(l == i || l == j || l == k) continue;
                            // cout << eps - distance(inter2, X[l]) << endl; 
                            if(distance(inter2, availableU[l]) < eps) flag_in = true; 
                        }
                        if(!flag_in){
                            success = false;
                            // std::cout << "パターン2" << std::endl;
                        }
                    }
                }   // kのループ

                if(!flag_cross_j){
                    bool flag_in = false;
                    quaternion inter = cal_inter_2ball(availableU[i], availableU[j], targetU, eps);
                    if(distance(inter, targetU) > eps) continue;   // ε-ballに含まれない
                    flag_cross_i = true;  
    #pragma omp parallel for
                    for(int k = 0; k < N; k++){
                        if(k == i || k == j) continue;
                        if(distance(inter, availableU[k]) < eps) flag_in = true;
                    }
                    if(!flag_in){
                        success = false;
                        // std::cout << "パターン3" << std::endl;
                    } 
                }
            }   // jのループ
            if(!flag_cross_i){
                quaternion inter = cal_inter_2ball(availableU[i], targetU, eps);
                if(!inter.is_unitary()){
                    continue;
                    success = false; 
                    // std::cout << "パターン4" << std::endl;
                }   // 一応書いてるけど、必ず交点は持つはず(持たないかも)

                bool flag_in = false;
                for(int j = 0; j < N; j++){
                    if(j == i) continue;
                    if(distance(inter, availableU[j]) < eps) flag_in = true;
                }
                if(!flag_in){
                    success = false; 
                    // std::cout << "パターン5" << std::endl;
                }
                
            }
        }   // iのループ
        
        if(success) return true;
        else return false;
    }
}




