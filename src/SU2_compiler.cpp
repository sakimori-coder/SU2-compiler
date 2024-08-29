#pragma once

#include <vector>
#include <string>
// #include <Eigen/Core>
#include "type.hpp"
#include "quaternion.hpp"
#include "rings.hpp"
#include "U2_ZOmega.hpp"
#include "grid_solver.hpp"
#include "grid_operator.hpp"
#include "ExactSynthesis.hpp"


namespace SU2_Compiler{

// ITYPE ITYPE_sqrt(ITYPE n){
//     if(n < 0) return -1;
//     ITYPE ok = n;
//     ITYPE ng = -1;
//     while(ok - ng > 1){
//         ITYPE mid = (ok + ng) / 2;
//         if(mid * mid >= n) ok = mid;
//         else ng = mid;
//     }
//     if(ok * ok == n) return ok;
//     else return -1;
// }


ITYPE ITYPE_sqrt(ITYPE n){
    double n_double = n;
    ITYPE n_sqrt = (ITYPE)(std::sqrt(n_double) + 0.5);
    if(n_sqrt * n_sqrt == n) return n_sqrt;
    else return -1;
}


ZRoot2 ZRoot2_sqrt(ZRoot2 n){
    ITYPE alpha = n.a;
    ITYPE beta = n.b;
    if(beta % 2 != 0) return {-1, 0}; 
    
    if(beta == 0){
        ITYPE a = ITYPE_sqrt(alpha);
        if(a != -1) return {a, 0};
        
        if(alpha % 2 == 0){
            ITYPE b = ITYPE_sqrt(alpha / 2);
            return {a, 0};
        }
    }else{
        ITYPE a = -1, b = -1;
        ITYPE tmp = ITYPE_sqrt(alpha * alpha - 2 * beta * beta);
        if(tmp != -1){
            if((alpha + tmp) % 2 == 0){
                ITYPE a_square1 = (alpha + tmp) / 2;
                ITYPE a_square2 = (alpha - tmp) / 2;
                ITYPE a1 = ITYPE_sqrt(a_square1);
                ITYPE a2 = ITYPE_sqrt(a_square2);
                if(a1 != -1) a = a1;
                if(a2 != -1) a = a2;
                if(a != -1){
                    if(beta % (2*a) == 0) b = beta / (2*a);
                }
            }
            if(a != -1 && b != -1){
                return {a, b};
            }
        }
    }

    return {-1, 0};
}


bool divisible(ZRoot2 x, ZRoot2 y){
    ZRoot2 numerator = x * conj(y);
    ITYPE denominator = y.norm();
    if(numerator.a % denominator == 0 && numerator.b % denominator == 0) return true;
    return false;
}


std::vector<ZRoot2> quadratic_equation_ZRoot2(ZRoot2 a, ZRoot2 b, ZRoot2 c){
    std::vector<ZRoot2> ret;
    ZRoot2 D_square = b * b - 4 * a * c;
    ZRoot2 D = ZRoot2_sqrt(D_square);
    if(D != ZRoot2(-1, 0)){
        ZRoot2 denominator = 2 * a;
        ZRoot2 numerator1 = -b + D;
        ZRoot2 numerator2 = -b - D;
        if(divisible(numerator1, denominator)){
            ret.push_back(numerator1 / denominator);
            ZRoot2 x = numerator1 / denominator;
            // std::cout << "左辺 " << a*x*x + b*x + c << std::endl;
        } 
        if(divisible(numerator2, denominator)){
            ret.push_back(numerator2 / denominator);
            ZRoot2 x = numerator2 / denominator;
            // std::cout << "左辺 " << a*x*x + b*x + c << std::endl;
        }
    }
    return ret;
}

std::string SU2_compiler(quaternion U, FTYPE eps){
    Vector4f p;
    p << 1, 0, 0, 0;
    FTYPE a = ((eps*eps) / 2);
    FTYPE b = sqrt(1 - (1 - eps*eps / 2) * (1 - eps*eps / 2));
    Matrix4f D, Delta, P;
    D << a*a, 0.0, 0.0, 0.0,
         0.0, b*b, 0.0, 0.0,
         0.0, 0.0, b*b, 0.0,
         0.0, 0.0, 0.0, b*b;
    std::cout << "D = \n" << D << std::endl;

    P << U.a,-U.b,-U.c,-U.d,
         U.b, U.a, U.d,-U.c,
         U.c,-U.d, U.a, U.b,
         U.d, U.c,-U.b, U.a;
    D = P * D * P.transpose();
    Delta = Matrix4f::Identity();
    p = P * p;

    Matrix4f X,Y;
    X << a, 0.0, 0.0, 0.0,
         0.0, b, 0.0, 0.0,
         0.0, 0.0, b, 0.0,
         0.0, 0.0, 0.0, b;
    X = P * X * P.transpose();
    X /= sqrt(sqrt(a*b*b*b));
    Y = Matrix4f::Identity();

    GridOp<4> G_T = LLL_ZRoot2<4>(X, Y, (FTYPE)0.99);
    GridOp<4> G = G_T.transpose();
    // G = GridOp<4>::Identity();
    GridOp<4> G_dot = conj(G);
    GridOp<4> G_inv = inverse_Special_GridOp(G);
    
    Matrix4f G_FTYPE = ZRoot2_to_FTYPE(G);
    Matrix4f G_dot_FTYPE = ZRoot2_to_FTYPE(G_dot);
    Matrix4f G_inv_FTYPE = ZRoot2_to_FTYPE(G_inv);

    Matrix4f GGT = G_FTYPE * G_FTYPE.transpose();
    std::cout << GGT << std::endl;


    // std::cout << determinant(G_inv) << std::endl;
    // std::cout << determinant(G_inv_dot) << std::endl;

    std::cout << "D = \n" << D << std::endl;
    p = G_FTYPE * p;
    D = G_FTYPE * D * G_FTYPE.transpose();
    Delta = G_dot_FTYPE * Delta * G_dot_FTYPE.transpose();

    std::cout << "D = \n" << D << std::endl;
    std::cout << "Δ = \n" << Delta << std::endl;

    std::array<std::pair<FTYPE, FTYPE>, 4> A, B;
    for(int i = 0; i < 4; i++){
        A[i] = {p(i) - sqrt(D(i,i)), p(i) + sqrt(D(i,i))};
        B[i] = {-sqrt(Delta(i,i)), sqrt(Delta(i,i))};
    }

    for(int k = 0; k < 100; k++){
        FTYPE sqrt2k = pow(sqrt2, (FTYPE)k);
        
        std::vector<ZRoot2> X = one_dim_grid_problem(A[0].first * sqrt2k, A[0].second * sqrt2k, 
                                                     B[0].first * sqrt2k, B[0].second * sqrt2k);

        std::vector<ZRoot2> Y = one_dim_grid_problem(A[1].first * sqrt2k, A[1].second * sqrt2k, 
                                                     B[1].first * sqrt2k, B[1].second * sqrt2k);

        std::vector<ZRoot2> Z = one_dim_grid_problem(A[2].first * sqrt2k, A[2].second * sqrt2k, 
                                                     B[2].first * sqrt2k, B[2].second * sqrt2k);

        std::vector<ZRoot2> W = one_dim_grid_problem(A[3].first * sqrt2k, A[3].second * sqrt2k, 
                                                     B[3].first * sqrt2k, B[3].second * sqrt2k);
        
        std::cout << "k = " << k << std::endl;
        std::cout << X.size() << " " << Y.size() << " " << Z.size() << " " << W.size() << std::endl;
 
        GridOp<4> D2 = G_inv.transpose() * G_inv;
        ZRoot2 key = {(1LL<<k), 0};
        std::vector<Eigen::Matrix<ZRoot2, 4, 1>> V_list;
#pragma omp parallel for
        for(auto &x : X){
#pragma omp parallel for
            for(auto &y : Y){
#pragma omp parallel for
                for(auto &z : Z){
                    ZRoot2 a = D2(3, 3);
                    ZRoot2 b = 2*D2(0,3)*x + 2*D2(1,3)*y + 2*D2(2,3)*z;
                    ZRoot2 c = D2(0,0)*x*x + D2(1,1)*y*y + D2(2,2)*z*z
                             + 2*D2(0,1)*x*y + 2*D2(0,2)*x*z + 2*D2(1,2)*y*z
                             - key;
                    
                    auto W_cand = quadratic_equation_ZRoot2(a,b,c);
                    for(auto &w : W_cand){
                        Eigen::Matrix<ZRoot2, 4, 1> xyzw;
                        xyzw << x, y, z, w;
                        xyzw = G_inv * xyzw;
                        ZRoot2 abs_xyzw = xyzw.dot(xyzw);
                        std::cout << "||v||^2 = " << abs_xyzw << std::endl;
                        ZRoot2 tmp = a * w*w + b * w + c;
#pragma parallel critical
                        V_list.push_back(xyzw);
                    }

                    // for(auto &w : W){
                    //     Eigen::Matrix<ZRoot2, 4, 1> xyzw;
                    //     xyzw << x, y, z, w;
                    //     xyzw = G_inv * xyzw;
                    //     ZRoot2 abs_xyzw = xyzw.dot(xyzw);
                    //     if(abs_xyzw == key){
                    //         ZRoot2 tmp = a * w*w + b * w + c;
                    //         std::cout << "2次方程式の左辺 " << tmp << std::endl;
                    //         return "success!";
                    //     }
                    // }
                }
            }
        }
        
        if(!V_list.empty()) return "success!";
    }
    
    return "failure";
}

}