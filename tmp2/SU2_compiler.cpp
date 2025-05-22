#include "su2compiler.hpp"

#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/LU>
#include "type.hpp"
#include "quaternion.hpp"
#include "rings.hpp"
#include "U2_ZOmega.hpp"
// #include "grid_solver.hpp"
// #include "grid_operator.hpp"
// #include "Unimodular.hpp"
#include "Unimodular_Z.hpp"
#include "ExactSynthesis.hpp"
#include <chrono>

using std::cout;
using std::endl;

namespace su2compiler{

ITYPE ITYPE_sqrt(ITYPE n){
    if(n < 0) return -1;
    ITYPE ok = n;
    ITYPE ng = -1;
    while(ok - ng > 1){
        ITYPE mid = (ok + ng) / 2;
        if(mid * mid >= n) ok = mid;
        else ng = mid;
    }
    if(ok * ok == n) return ok;
    else return -1;
}


// ITYPE ITYPE_sqrt(ITYPE n){
//     double n_double = n;
//     ITYPE n_sqrt = (ITYPE)(std::sqrt(n_double) + 0.5);
//     if(n_sqrt * n_sqrt == n) return n_sqrt;
//     else return -1;
// }


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





std::chrono::system_clock::time_point start, end;
double time;



// std::string su2compiler(quaternion U, FTYPE eps, int l){
//     Matrix4f P;
//     P << U.a,-U.b,-U.c,-U.d,
//          U.b, U.a, U.d,-U.c,
//          U.c,-U.d, U.a, U.b,
//          U.d, U.c,-U.b, U.a;

//     for(int k = 0; k < 100; k++){
//         // std::cout << "k = " << k << std::endl;
//         FTYPE r = pow(sqrt2, (FTYPE)k);
//         std::vector<std::string> V_list;


//         start = std::chrono::system_clock::now();
//         // std::cout << "start" << std::endl;
//         // 楕円Aを定義
//         Matrix4f D;
//         Vector4f p_list[2];
//         FTYPE a = r * (1 - sqrt(1 - eps*eps));
//         FTYPE b = r * eps;
//         D << a*a, 0.0, 0.0, 0.0,
//             0.0, b*b, 0.0, 0.0,
//             0.0, 0.0, b*b, 0.0,
//             0.0, 0.0, 0.0, b*b;
//         D = P * D * P.transpose();
//         p_list[0] << r * sqrt(1 - eps*eps), 0, 0, 0;
//         p_list[0] = P * p_list[0];
//         for(int j = 0; j < 4; j++) p_list[1](j) = p_list[0](j) - inv_sqrt2;

//         // 基底行列の定義
//         Matrix4f X;
//         X << a, 0.0, 0.0, 0.0,
//                 0.0, b, 0.0, 0.0,
//                 0.0, 0.0, b, 0.0,
//                 0.0, 0.0, 0.0, b;
//         X = P * X * P.transpose();


//         // 楕円Bを定義
//         Matrix4f Delta, Q;
//         Vector4f q_list[2];
//         Delta = r*r * Matrix4f::Identity();
//         q_list[0] << 0.0, 0.0, 0.0, 0.0;
//         for(int j = 0; j < 4; j++) q_list[1](j) = q_list[0](j) + inv_sqrt2;

//         // 基底行列の定義
//         Matrix4f Y;
//         Y = r * Matrix4f::Identity();

//         // ユニモジュラ行列の計算と作用
//         start = std::chrono::system_clock::now();
//         UnimodularMatrix<4> G_T = LLL_ZRoot2<4>(X, Y, (FTYPE)0.99);
//         // G_T = UnimodularMatrix<4>::Identity();
//         UnimodularMatrix<4> G = G_T.transpose();
//         UnimodularMatrix<4> G_dot = conj(G);
//         UnimodularMatrix<4> G_inv = inverse_UnimodularMatrix(G);
        
//         Matrix4f G_FTYPE = ZRoot2_to_FTYPE(G);
//         Matrix4f G_dot_FTYPE = ZRoot2_to_FTYPE(G_dot);
//         Matrix4f G_inv_FTYPE = ZRoot2_to_FTYPE(G_inv);

//         for(int i = 0; i < 2; i++){
//             p_list[i] = G_FTYPE * p_list[i];
//             q_list[i] = G_dot_FTYPE * q_list[i];
//         }
//         D = G_FTYPE * D * G_FTYPE.transpose();
//         Delta = G_dot_FTYPE * Delta * G_dot_FTYPE.transpose();

//         end = std::chrono::system_clock::now();
//         time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
//         // std::cout << "time1 " << time << "[ms]" << std::endl;
        
        
//         start = std::chrono::system_clock::now();
//         for(int u_omega = 0; u_omega < 2; u_omega++){
//             for(int t_omega = 0; t_omega < 2; t_omega++){
//                 Vector4f p = p_list[u_omega];
//                 Vector4f q = q_list[t_omega];
//                 std::array<std::pair<FTYPE, FTYPE>, 4> A, B;
//                 for(int j = 0; j < 4; j++){
//                     A[j] = {p(j) - sqrt(D(j,j)), p(j) + sqrt(D(j,j))};
//                     B[j] = {q(j) - sqrt(Delta(j,j)), q(j) + sqrt(Delta(j,j))};
//                     // std::cout << A[j].second - A[j].first << std::endl;
//                     // std::cout << B[j].second - B[j].first << std::endl;
//                     // std::cout << std::endl;
//                 }


//                 std::vector<ZRoot2> X = one_dim_grid_problem(A[0].first, A[0].second, 
//                                                             B[0].first, B[0].second);

//                 std::vector<ZRoot2> Y = one_dim_grid_problem(A[1].first, A[1].second, 
//                                                             B[1].first, B[1].second);

//                 std::vector<ZRoot2> Z = one_dim_grid_problem(A[2].first, A[2].second, 
//                                                             B[2].first, B[2].second);

//                 std::vector<ZRoot2> W = one_dim_grid_problem(A[3].first, A[3].second, 
//                                                             B[3].first, B[3].second);


//                 // std::cout << X.size() << " " << Y.size() << " " << Z.size() << " " << W.size() << std::endl;
//                 // std::cout << "解の個数 : " << X.size() * Y.size() * Z.size() * W.size() << std::endl;
//                 // if(X.size() * Y.size() * Z.size() * W.size() > 100000){
//                 //     std::cout << "解の個数 : " << X.size() * Y.size() * Z.size() * W.size() << std::endl;
//                 // }

//                 // i=0の場合は(√2x)^2 + (√2y)^2 + (√2z)^2 + (√2w)^2 = 2^{k+1}
//                 // i=1の場合は(√2x + 1)^2 + (√2y + 1)^2 + (√2z + 1)^2 + (√2w + 1)^2 = 2^{k+1}
//                 // を満たす(x,y,z,w)を探す。
//                 for(auto &x : X) x *= ZRoot2(0, 1);
//                 for(auto &y : Y) y *= ZRoot2(0, 1);
//                 for(auto &z : Z) z *= ZRoot2(0, 1);
//                 for(auto &w : W) w *= ZRoot2(0, 1);

//                 Eigen::Vector<ZRoot2, 4> shift;
//                 if(u_omega == 0){
//                     shift(0) = ZRoot2(0,0);
//                     shift(1) = ZRoot2(0,0);
//                 }
//                 else{
//                     shift(0) = ZRoot2(1,0);
//                     shift(1) = ZRoot2(1,0);
//                 }
                
//                 if(t_omega == 0){
//                     shift(2) = ZRoot2(0,0);
//                     shift(3) = ZRoot2(0,0);
//                 }
//                 else{
//                     shift(2) = ZRoot2(1,0);
//                     shift(3) = ZRoot2(1,0);
//                 }
//                 shift = G * shift;
//                 for(auto &x : X) x += shift(0);
//                 for(auto &y : Y) y += shift(1);
//                 for(auto &z : Z) z += shift(2);
//                 for(auto &w : W) w += shift(3);

                
//                 ZRoot2 key = {((ITYPE)1<<(k+1)), 0};
//                 UnimodularMatrix<4> H;
//                 H = G_inv.transpose() * G_inv;
// #pragma omp parallel for collapse(3)
//                 for(auto &x : X){
//                     for(auto &y : Y){
//                         for(auto &z : Z){
//                             ZRoot2 a = H(3, 3);
//                             ZRoot2 b = 2*H(0,3)*x + 2*H(1,3)*y + 2*H(2,3)*z;
//                             ZRoot2 c = H(0,0)*x*x + H(1,1)*y*y + H(2,2)*z*z
//                                        + 2*H(0,1)*x*y + 2*H(0,2)*x*z + 2*H(1,2)*y*z
//                                        - key;
//                             auto W_cand = quadratic_equation_ZRoot2(a,b,c);
//                             // std::cout << W_cand.size() << std::endl;
//                             for(auto &w : W_cand){
//                                 Eigen::Matrix<ZRoot2, 4, 1> xyzw;
//                                 xyzw << x, y, z, w;
//                                 xyzw = G_inv * xyzw;
//                                 // √2倍されてるから
//                                 ITYPE a0 = xyzw(0).b;
//                                 ITYPE b0 = (xyzw(0).a - (ITYPE)u_omega) / 2;
//                                 ITYPE a1 = xyzw(1).b;
//                                 ITYPE b1 = (xyzw(1).a - (ITYPE)u_omega) / 2;
//                                 ITYPE a2 = xyzw(2).b;
//                                 ITYPE b2 = (xyzw(2).a - (ITYPE)t_omega) / 2;
//                                 ITYPE a3 = xyzw(3).b;
//                                 ITYPE b3 = (xyzw(3).a - (ITYPE)t_omega) / 2;

//                                 ZOmega u,t;
//                                 if(u_omega == 0) u = {-b0+b1, a1, b0+b1, a0};
//                                 else u = {-b0+b1, a1, b0+b1+1, a0};
                                    
//                                 if(t_omega == 0) t = {-b2+b3, a3, b2+b3, a2};
//                                 else t = {-b2+b3, a3, b2+b3+1, a2};
    
//                                 if(u.norm() + t.norm() != ZRoot2((ITYPE)1 << k, 0)) continue;

//                                 U2_ZOmega V(u, t, l, k);
//                                 // std::cout << "dist " << distance(U, to_quaternion(V)) << std::endl;
//                                 if(distance(U, to_quaternion(V)) < eps){
//                                     std::string V_str = ExactSynthesis(V);
// #pragma omp critical
//                                     if(V_str == "SHTSHTSHTHTSHTSHTHTHTHTHTSHTSHTSHTHTSHTSHTSHTHTSHTSHTSHTHTHTSHTSHTHTHTHTHTSHTHTHTSHTHTSHTHTHTSHTHTSHTHTHTSHTHTSHTSHTHTHTHTSHTSHTHTHTSHTSHTSHTSSHSS"){
//                                         std::cout << k << std::endl;
//                                         std::cout << a0 << " " << b0 << std::endl;
//                                         std::cout << a1 << " " << b1 << std::endl;
//                                         std::cout << a2 << " " << b2 << std::endl;
//                                         std::cout << a3 << " " << b3 << std::endl;
//                                     }
//                                     V_list.push_back(V_str);
//                                 }
//                             } // W
//                         } // Z
//                     } // Y 
//                 } // X
//             } // t_omega
//         } // u_omega
//         end = std::chrono::system_clock::now();
//         time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
//         // std::cout << "time2 " << time << "[ms]" << std::endl;
        
//         if(!V_list.empty()){
//             auto comp = [](const std::string& x, const std::string &y){
//                 return std::count(x.begin(), x.end(), 'T') < std::count(y.begin(), y.end(), 'T');
//             };
//             return *min_element(V_list.begin(), V_list.end(), comp);
//         }
//     }
//     return "Failure";
// }



// std::vector<U2_ZOmega> enum_u_t(quaternion U, FTYPE eps, int k_max, int l){
//     quaternion phase(1,0,0,0);
//     quaternion sqrt_omega_FTYPE(sqrt_omega.real(), sqrt_omega.imag(), 0, 0);
//     for(int i = 0; i < l; i++) phase *= sqrt_omega_FTYPE;
//     quaternion U_prime = U * phase;

//     Matrix4f P;
//     P << U_prime.a,-U_prime.b,-U_prime.c,-U_prime.d,
//          U_prime.b, U_prime.a, U_prime.d,-U_prime.c,
//          U_prime.c,-U_prime.d, U_prime.a, U_prime.b,
//          U_prime.d, U_prime.c,-U_prime.b, U_prime.a;

//     FTYPE r = pow(sqrt2, (FTYPE)k_max);
//     std::vector<U2_ZOmega> V_list;


//     start = std::chrono::system_clock::now();
//     // std::cout << "start" << std::endl;
//     // 楕円Aを定義
//     Matrix4f D;
//     Vector4f p_list[2];
//     FTYPE a = r * (1 - sqrt(1 - eps*eps));
//     FTYPE b = r * eps;
//     D << a*a, 0.0, 0.0, 0.0,
//         0.0, b*b, 0.0, 0.0,
//         0.0, 0.0, b*b, 0.0,
//         0.0, 0.0, 0.0, b*b;
//     D = P * D * P.transpose();
//     p_list[0] << r * sqrt(1 - eps*eps), 0, 0, 0;
//     p_list[0] = P * p_list[0];
//     for(int j = 0; j < 4; j++) p_list[1](j) = p_list[0](j) - inv_sqrt2;

//     // 基底行列の定義
//     Matrix4f X;
//     X << a, 0.0, 0.0, 0.0,
//             0.0, b, 0.0, 0.0,
//             0.0, 0.0, b, 0.0,
//             0.0, 0.0, 0.0, b;
//     X = P * X * P.transpose();


//     // 楕円Bを定義
//     Matrix4f Delta, Q;
//     Vector4f q_list[2];
//     Delta = r*r * Matrix4f::Identity();
//     q_list[0] << 0.0, 0.0, 0.0, 0.0;
//     for(int j = 0; j < 4; j++) q_list[1](j) = q_list[0](j) + inv_sqrt2;

//     // 基底行列の定義
//     Matrix4f Y;
//     Y = r * Matrix4f::Identity();

//     // ユニモジュラ行列の計算と作用
//     start = std::chrono::system_clock::now();
//     UnimodularMatrix<4> G_T = LLL_ZRoot2<4>(X, Y, (FTYPE)0.99);
//     // G_T = UnimodularMatrix<4>::Identity();
//     UnimodularMatrix<4> G = G_T.transpose();
//     UnimodularMatrix<4> G_dot = conj(G);
//     UnimodularMatrix<4> G_inv = inverse_UnimodularMatrix(G);
    
//     Matrix4f G_FTYPE = ZRoot2_to_FTYPE(G);
//     Matrix4f G_dot_FTYPE = ZRoot2_to_FTYPE(G_dot);
//     Matrix4f G_inv_FTYPE = ZRoot2_to_FTYPE(G_inv);

//     for(int i = 0; i < 2; i++){
//         p_list[i] = G_FTYPE * p_list[i];
//         q_list[i] = G_dot_FTYPE * q_list[i];
//     }
//     D = G_FTYPE * D * G_FTYPE.transpose();
//     Delta = G_dot_FTYPE * Delta * G_dot_FTYPE.transpose();

//     end = std::chrono::system_clock::now();
//     time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
//     // std::cout << "time1 " << time << "[ms]" << std::endl;
    
    
//     start = std::chrono::system_clock::now();
//     for(int u_omega = 0; u_omega < 2; u_omega++){
//         for(int t_omega = 0; t_omega < 2; t_omega++){
//             Vector4f p = p_list[u_omega];
//             Vector4f q = q_list[t_omega];
//             std::array<std::pair<FTYPE, FTYPE>, 4> A, B;
//             for(int j = 0; j < 4; j++){
//                 A[j] = {p(j) - sqrt(D(j,j)), p(j) + sqrt(D(j,j))};
//                 B[j] = {q(j) - sqrt(Delta(j,j)), q(j) + sqrt(Delta(j,j))};
//                 // std::cout << A[j].second - A[j].first << std::endl;
//                 // std::cout << B[j].second - B[j].first << std::endl;
//                 // std::cout << std::endl;
//             }


//             std::vector<ZRoot2> X = one_dim_grid_problem(A[0].first, A[0].second, 
//                                                         B[0].first, B[0].second);

//             std::vector<ZRoot2> Y = one_dim_grid_problem(A[1].first, A[1].second, 
//                                                         B[1].first, B[1].second);

//             std::vector<ZRoot2> Z = one_dim_grid_problem(A[2].first, A[2].second, 
//                                                         B[2].first, B[2].second);

//             std::vector<ZRoot2> W = one_dim_grid_problem(A[3].first, A[3].second, 
//                                                         B[3].first, B[3].second);


//             // std::cout << X.size() << " " << Y.size() << " " << Z.size() << " " << W.size() << std::endl;
//             // std::cout << "解の個数 : " << X.size() * Y.size() * Z.size() * W.size() << std::endl;
//             // if(X.size() * Y.size() * Z.size() * W.size() > 100000){
//             //     std::cout << "解の個数 : " << X.size() * Y.size() * Z.size() * W.size() << std::endl;
//             // }

//             // i=0の場合は(√2x)^2 + (√2y)^2 + (√2z)^2 + (√2w)^2 = 2^{k+1}
//             // i=1の場合は(√2x + 1)^2 + (√2y + 1)^2 + (√2z + 1)^2 + (√2w + 1)^2 = 2^{k+1}
//             // を満たす(x,y,z,w)を探す。
//             for(auto &x : X) x *= ZRoot2(0, 1);
//             for(auto &y : Y) y *= ZRoot2(0, 1);
//             for(auto &z : Z) z *= ZRoot2(0, 1);
//             for(auto &w : W) w *= ZRoot2(0, 1);

//             Eigen::Vector<ZRoot2, 4> shift;
//             if(u_omega == 0){
//                 shift(0) = ZRoot2(0,0);
//                 shift(1) = ZRoot2(0,0);
//             }
//             else{
//                 shift(0) = ZRoot2(1,0);
//                 shift(1) = ZRoot2(1,0);
//             }
            
//             if(t_omega == 0){
//                 shift(2) = ZRoot2(0,0);
//                 shift(3) = ZRoot2(0,0);
//             }
//             else{
//                 shift(2) = ZRoot2(1,0);
//                 shift(3) = ZRoot2(1,0);
//             }
//             shift = G * shift;
//             for(auto &x : X) x += shift(0);
//             for(auto &y : Y) y += shift(1);
//             for(auto &z : Z) z += shift(2);
//             for(auto &w : W) w += shift(3);

            
//             ZRoot2 key = {((ITYPE)1<<(k_max+1)), 0};
//             UnimodularMatrix<4> H;
//             H = G_inv.transpose() * G_inv;
// #pragma omp parallel for collapse(3)
//             for(auto &x : X){
//                 for(auto &y : Y){
//                     for(auto &z : Z){
//                         ZRoot2 a = H(3, 3);
//                         ZRoot2 b = 2*H(0,3)*x + 2*H(1,3)*y + 2*H(2,3)*z;
//                         ZRoot2 c = H(0,0)*x*x + H(1,1)*y*y + H(2,2)*z*z
//                                     + 2*H(0,1)*x*y + 2*H(0,2)*x*z + 2*H(1,2)*y*z
//                                     - key;
//                         auto W_cand = quadratic_equation_ZRoot2(a,b,c);
//                         // std::cout << W_cand.size() << std::endl;
//                         for(auto &w : W_cand){
//                             Eigen::Matrix<ZRoot2, 4, 1> xyzw;
//                             xyzw << x, y, z, w;
//                             xyzw = G_inv * xyzw;
//                             // √2倍されてるから
//                             ITYPE a0 = xyzw(0).b;
//                             ITYPE b0 = (xyzw(0).a - (ITYPE)u_omega) / 2;
//                             ITYPE a1 = xyzw(1).b;
//                             ITYPE b1 = (xyzw(1).a - (ITYPE)u_omega) / 2;
//                             ITYPE a2 = xyzw(2).b;
//                             ITYPE b2 = (xyzw(2).a - (ITYPE)t_omega) / 2;
//                             ITYPE a3 = xyzw(3).b;
//                             ITYPE b3 = (xyzw(3).a - (ITYPE)t_omega) / 2;

//                             ZOmega u,t;
//                             if(u_omega == 0) u = {-b0+b1, a1, b0+b1, a0};
//                             else u = {-b0+b1, a1, b0+b1+1, a0};
                                
//                             if(t_omega == 0) t = {-b2+b3, a3, b2+b3, a2};
//                             else t = {-b2+b3, a3, b2+b3+1, a2};

//                             if(u.norm() + t.norm() != ZRoot2((ITYPE)1 << k_max, 0)) continue;

//                             U2_ZOmega V(u, t, l, k_max);
//                             // std::cout << "dist " << distance(U, to_quaternion(V)) << std::endl;
//                             if(distance(U, to_quaternion(V)) < eps){
// #pragma omp critical
//                                 V_list.push_back(V);
//                             }
//                         } // W
//                     } // Z
//                 } // Y 
//             } // X
//         } // t_omega
//     } // u_omega
//     end = std::chrono::system_clock::now();
//     time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
//     // std::cout << "time2 " << time << "[ms]" << std::endl;
    
//     return V_list;
// }



// template<int M>
// std::vector<Eigen::Vector<ITYPE, -1>> Enumerate_Integer_Points(
//     Eigen::Matrix<FTYPE, M, M> Q, 
//     Eigen::Vector<FTYPE, M> p)
// {
//     std::vector<Eigen::Matrix<FTYPE, -1, -1>> Q_sub(M+1), Q_sub_inv(M+1);
//     for(int i = 1; i <= M; i++){
//         Q_sub[i] = Q.block(M-i, M-i, i, i);
//         Q_sub_inv[i] = Q.block(M-i, M-i, i, i).inverse();
//     }
        
//     // (x-p)^T Q (x-p) <= cを満たす整数ベクトルxを列挙
//     std::function<std::vector<VectorXi>(int, VectorXf, FTYPE)> Enum;
//     Enum = [&Enum, &Q_sub, &Q_sub_inv](int N, VectorXf p, FTYPE c){
//         std::vector<VectorXi> ret;
//         if(c < 0) return ret;

//         MatrixXf Q = Q_sub[N];
//         if(N == 1){
//             ITYPE x0 = (ITYPE)ceil(p(0) - sqrt(Q_sub_inv[N](0,0) * c));
//             ITYPE x1 = (ITYPE)floor(p(0) + sqrt(Q_sub_inv[N](0,0) * c));
//             for(ITYPE x = x0; x <= x1; x++){
//                 VectorXi sol(N);
//                 sol(0) = x;
//                 ret.push_back(sol);
//             }
//             return ret;
//         }
        
//         ITYPE x0 = (ITYPE)ceil(p(0) - sqrt(Q_sub_inv[N](0,0) * c));
//         ITYPE x1 = (ITYPE)floor(p(0) + sqrt(Q_sub_inv[N](0,0) * c));
//         for(ITYPE x = x0; x <= x1; x++){
//             VectorXi sol(N);
//             sol(0) = x;
            
//             MatrixXf A_inv = Q_sub_inv[N-1];
//             VectorXf b = (x - p(0)) * Q.col(0).tail(N-1);
//             FTYPE next_c = c;
//             VectorXf next_p = p.tail(N-1);
//             next_p -= A_inv * b;
//             next_c -= Q(0,0) * (x - p(0)) * (x - p(0));
//             next_c += b.dot(A_inv * b);
//             // std::cout << N << " " << c << std::endl;
//             auto X_tail = Enum(N-1, next_p, next_c);
//             for(auto x_tail : X_tail){
//                 sol.tail(N-1) = x_tail;
//                 ret.push_back(sol);
//             }
//         } 
//         return ret;
//     };

//     auto ret = Enum(M, p, 1.0);
//     return ret;
// }


// (x-p)^T D^{-1} (x-p) <= cを満たす整数点xを列挙
std::vector<VectorXi> Enumerate_Integer_Points(
    MatrixXf Q, 
    VectorXf p,
    FTYPE c)
{
    MatrixXf Q_inv = Q.inverse();
    std::vector<VectorXi> ret;
    if(N == 1){
        ITYPE x0 = (ITYPE)ceil(p(0) - sqrt(Q_inv(0,0) * c));
        ITYPE x1 = (ITYPE)floor(p(0) + sqrt(Q_inv(0,0) * c));
        for(ITYPE x = x0; x <= x1; x++){
            VectorXi sol;
            sol(0) = x;
            ret.push_back(sol);
        }
    }else{
        Eigen::LLT<MatrixXf> llt_of_Qinv(Q_inv);
        MatrixXf B = llt_of_Qinv.matrixL().transpose();
        UnimodularMatrix<N> U_T = LLL(B, 0.5);
        UnimodularMatrix<N> U = U_T.transpose();
        UnimodularMatrix<N> U_inv = inverse_UnimodularMatrix(U);

        Eigen::Matrix<FTYPE, N, N> U_FTYPE;
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                U_FTYPE(i,j) = (FTYPE)U(i,j);
            }
        }

        p = U_FTYPE * p;
        Q_inv = U_FTYPE * Q_inv * U_FTYPE.transpose();
        
        ITYPE x0 = (ITYPE)ceil(p(0) - sqrt(Q_inv(0,0) * c));
        ITYPE x1 = (ITYPE)floor(p(0) + sqrt(Q_inv(0,0) * c));
        for(ITYPE x = x0; x <= x1; x++){
            Eigen::Vector<ITYPE, N> sol;
            sol(0) = x;
            
            Eigen::Matrix<FTYPE, N-1, N-1> A_inv = Q.block(1,1, N-1,N-1);
            Eigen::Vector<FTYPE, N-1> b = (x - p(0)) * Q.col(0).tail(N-1);
            FTYPE next_c = c;
            Eigen::Vector<FTYPE, N-1> next_p = p.tail(N-1);
            next_p -= A_inv * b;
            next_c -= Q(0,0) * (x - p(0)) * (x - p(0));
            next_c += b.dot(A_inv * b);
            // std::cout << N << " " << c << std::endl;
            auto X_tail = Enumerate_Integer_Points<N-1>(Q.block(1,1, N-1,N-1), next_p, next_c);
            for(auto x_tail : X_tail){
                sol.tail(N-1) = x_tail;
                sol = U_inv * sol;
                ret.push_back(sol);
            }
        }
    }

    return ret;
}


std::string su2compiler(quaternion U, FTYPE eps, int l){

    // 楕円の固有ベクトル
    Matrix8f P;
    P << U.a,-U.b,-U.c,-U.d, 0, 0, 0, 0,
         U.b, U.a, U.d,-U.c, 0, 0, 0, 0,
         U.c,-U.d, U.a, U.b, 0, 0, 0, 0,
         U.d, U.c,-U.b, U.a, 0, 0, 0, 0,
         0, 0, 0, 0,         1, 0, 0, 0,
         0, 0, 0, 0,         0, 1, 0, 0,
         0, 0, 0, 0,         0, 0, 1, 0,
         0, 0, 0, 0,         0, 0, 0, 1;         

    Matrix8f Sigma, Sigma_inv;
    Sigma << 1, 1/sqrt2, 0,-1/sqrt2, 0, 0, 0, 0,
             0, 1/sqrt2, 1, 1/sqrt2, 0, 0, 0, 0,
             0, 0, 0, 0, 1, 1/sqrt2, 0,-1/sqrt2,
             0, 0, 0, 0, 0, 1/sqrt2, 1, 1/sqrt2,
             1,-1/sqrt2, 0, 1/sqrt2, 0, 0, 0, 0,
             0,-1/sqrt2, 1,-1/sqrt2, 0, 0, 0, 0,
             0, 0, 0, 0, 1,-1/sqrt2, 0, 1/sqrt2,
             0, 0, 0, 0, 0,-1/sqrt2, 1,-1/sqrt2;

    Sigma_inv = Sigma.transpose();
    Sigma_inv /= 2.0;


    int k = 0;
    while(true){
        std::cout << "k = " << k << std::endl;
        FTYPE r = pow(sqrt2, (FTYPE)k);
        std::vector<std::string> V_list;


        start = std::chrono::system_clock::now();
        // std::cout << "start" << std::endl;
        // 楕円Aを定義
        Matrix8f B, D;
        Vector8f p;
        FTYPE e1 = r * (1 - sqrt(1 - eps*eps));
        FTYPE e2 = r * eps;
        B <<  e1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0,  e2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0,  e2, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0,  e2, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0,   r, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0,   r, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   r, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   r;
        B *= sqrt2;
        B = P * B * P.transpose();
        B = B * Sigma_inv.transpose();

        D = B.transpose() * B;
        p << r * sqrt(1 - eps*eps), 0, 0, 0, 0, 0, 0, 0;
        p = P * p;
        p = Sigma_inv * p;


        // ユニモジュラ行列の計算と作用
        start = std::chrono::system_clock::now();
        UnimodularMatrix<8> G_T = LLL<8>(B, 0.5);
        // G_T = UnimodularMatrix<8>::Identity();
        end = std::chrono::system_clock::now();
        time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        cout << "time : " << time << "[ms]" << endl;
        UnimodularMatrix<8> G = G_T.transpose();
        UnimodularMatrix<8> G_inv = inverse_UnimodularMatrix(G);

        Matrix8f G_FTYPE, G_inv_FTYPE;
        for(int i = 0; i < 8; i++){
            for(int j = 0; j < 8; j++){
                G_FTYPE(i,j) = (FTYPE)G(i,j);
                G_inv_FTYPE(i,j) = (FTYPE)G_inv(i,j);
            }
        }

        p = G_FTYPE * p;
        D = G_FTYPE * D * G_FTYPE.transpose();

        // std::cout << D << std::endl;
        // std::cout << (G.transpose() * G) << std::endl;
        // std::cout << (G_inv.transpose() * G_inv) << std::endl;
        // std::cout << Sigma * G_FTYPE.transpose() * Sigma_inv << std::endl;
        

        // end = std::chrono::system_clock::now();
        // time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
        // // std::cout << "time1 " << time << "[ms]" << std::endl;
        
        std::vector<ITYPE> Start(8), End(8);
        long long Total = 1;
        for(int i = 0; i < 8; i++){
            Start[i] = (ITYPE)ceil(p(i) - sqrt(D(i,i)));
            End[i] = (ITYPE)floor(p(i) + sqrt(D(i,i)));
            // std::cout << Start[i] << " " << End[i] << std::endl;
            // std::cout << End[i] - Start[i] << std::endl;
            Total *= (End[i] - Start[i] + 1);
        }
        std::cout << "Total=" << Total << std::endl;

        // std::cout << std::endl;
        // for(int i = 0; i < 8; i++){
        //     cout << i << ":" << Start[i] << " " << End[i] << endl;
        // }

        int cnt0 = 0, cnt1 = 0;

        Matrix8f D_inv = D.inverse();
        // auto solutions = Enum(8, D_inv, p, 1.0);
        auto solutions = Enumerate_Integer_Points<8>(D_inv, p, 1.0);
        std::cout << "End enumeration" << std::endl;
        for(auto v : solutions){
            // Vector8f v_FTYPE;
            // for(int i = 0; i < 8; i++) v_FTYPE(i) = (FTYPE)v(i);
            // FTYPE value = (v_FTYPE-p).dot(D_inv * (v_FTYPE-p));
            FTYPE value = 1.0;

            v = G_inv * v;
            ZOmega u(v(3), v(2), v(1), v(0));
            ZOmega t(v(7), v(6), v(5), v(4));

            U2_ZOmega V(u, t, l, k);
            // std::cout << "dist " << distance(U, to_quaternion(V)) << std::endl;
            ZRoot2 N = conj(u.norm()) + conj(t.norm());
            // N = u.norm() + t.norm();
            // std::cout << value << std::endl;
            // if(ZRoot2_to_FTYPE(N) > pow((FTYPE)2, k)){
            if(value > 1.0){
                cnt0++;
            }else{
                cnt1++;
            }

            if(u.norm() + t.norm() == ZRoot2((ITYPE)1 << k, 0)){
                std::cout << v.dot(v) << std::endl;
                std::cout << "dist " << distance(U, to_quaternion(V)) << std::endl;
                std::cout << "value " << value << std::endl;
                if(distance(U, to_quaternion(V)) < eps){
                    std::string V_str = ExactSynthesis(V);
                    V_list.push_back(V_str);
                }
            }
        }


        // for(ITYPE x0 = Start[0]; x0 <= End[0]; x0++){
        // for(ITYPE x1 = Start[1]; x1 <= End[1]; x1++){
        // for(ITYPE x2 = Start[2]; x2 <= End[2]; x2++){
        // for(ITYPE x3 = Start[3]; x3 <= End[3]; x3++){
        // for(ITYPE x4 = Start[4]; x4 <= End[4]; x4++){
        // for(ITYPE x5 = Start[5]; x5 <= End[5]; x5++){
        // for(ITYPE x6 = Start[6]; x6 <= End[6]; x6++){
        // for(ITYPE x7 = Start[7]; x7 <= End[7]; x7++){
        //     Vector8I v;
        //     v << x0, x1, x2, x3, x4, x5, x6, x7;
            
        //     // Vector8f v_FTYPE;
        //     // for(int i = 0; i < 8; i++) v_FTYPE(i) = (FTYPE)v(i);
        //     // FTYPE value = (v_FTYPE-p).dot(D_inv * (v_FTYPE-p));
        //     FTYPE value = 1.0;

        //     v = G_inv * v;
        //     ZOmega u(v(3), v(2), v(1), v(0));
        //     ZOmega t(v(7), v(6), v(5), v(4));

        //     U2_ZOmega V(u, t, l, k);
        //     // std::cout << "dist " << distance(U, to_quaternion(V)) << std::endl;
        //     ZRoot2 N = conj(u.norm()) + conj(t.norm());
        //     // N = u.norm() + t.norm();
        //     // std::cout << value << std::endl;
        //     // if(ZRoot2_to_FTYPE(N) > pow((FTYPE)2, k)){
        //     if(value > 1.0){
        //         cnt0++;
        //     }else{
        //         cnt1++;
        //     }



        //     if(u.norm() + t.norm() == ZRoot2((ITYPE)1 << k, 0)){
        //         std::cout << v.dot(v) << std::endl;
        //         std::cout << "dist " << distance(U, to_quaternion(V)) << std::endl;
        //         std::cout << "value " << value << std::endl;
        //         if(distance(U, to_quaternion(V)) < eps){
        //             std::string V_str = ExactSynthesis(V);
        //             V_list.push_back(V_str);
        //         }
        //     }
        // }
        // }
        // }
        // }
        // }
        // }
        // }
        // }
        
        std::cout << "cnt0:" << cnt0 << endl;
        std::cout << "cnt1:" << cnt1 << endl;

        if(!V_list.empty()){
            auto comp = [](const std::string& x, const std::string &y){
                return std::count(x.begin(), x.end(), 'T') < std::count(y.begin(), y.end(), 'T');
            };
            return *min_element(V_list.begin(), V_list.end(), comp);
        }
        k++;
        std::cout << "終了" << std::endl;
    }
    return "Failure";
}





// std::string su2compiler(quaternion U, FTYPE eps, int l){
//     Matrix4f P;
//     P << U.a,-U.b,-U.c,-U.d,
//          U.b, U.a, U.d,-U.c,
//          U.c,-U.d, U.a, U.b,
//          U.d, U.c,-U.b, U.a;

//     for(int k = 0; k < 100; k++){
//         std::cout << "k = " << k << std::endl;
//         FTYPE r = pow(sqrt2, (FTYPE)k);
//         std::vector<std::string> V_list;

//         FTYPE delta = min((FTYPE)1.0, 2.0*pow(r, -(FTYPE)1.6) / eps);
//         std::cout << "delta = " << delta << std::endl;
//         std::cout << "a = " << sqrt(8.0/(PI*PI) * (1 - pow(1-delta*delta, 1.0/(FTYPE)3.0))) << std::endl;
//         int N = (int)ceil(1.0/ sqrt(8.0/(PI*PI) * (1 - pow(1-delta*delta, 1.0/(FTYPE)3.0))));
//         std::cout << "N = " << N << std::endl;

//         std::vector<FTYPE> theta1_list(N), theta2_list(N), theta3_list(2*N);
//         for(int i = 0; i < N; i++) theta1_list[i] = PI * ((FTYPE)i / N);
//         for(int i = 0; i < N; i++) theta2_list[i] = PI * ((FTYPE)i / N);
//         for(int i = 0; i < 2*N; i++) theta3_list[i] = PI * ((FTYPE)i / N);

// #pragma omp parallel for collapse(3)
//         for(auto &theta1 : theta1_list){
//             for(auto &theta2 : theta2_list){
//                 for(auto &theta3 : theta3_list){
//                     start = std::chrono::system_clock::now();
//                     // std::cout << "start" << std::endl;
//                     // 楕円Aを定義
//                     Matrix4f D;
//                     Vector4f p_list[2];
//                     FTYPE a = r * (1 - sqrt(1 - eps*eps));
//                     FTYPE b = r * eps;
//                     D << a*a, 0.0, 0.0, 0.0,
//                         0.0, b*b, 0.0, 0.0,
//                         0.0, 0.0, b*b, 0.0,
//                         0.0, 0.0, 0.0, b*b;
//                     D = P * D * P.transpose();
//                     p_list[0] << r * sqrt(1 - eps*eps), 0, 0, 0;
//                     p_list[0] = P * p_list[0];
//                     for(int j = 0; j < 4; j++) p_list[1](j) = p_list[0](j) - inv_sqrt2;

//                     // 基底行列の定義
//                     Matrix4f X;
//                     X << a, 0.0, 0.0, 0.0,
//                          0.0, b, 0.0, 0.0,
//                          0.0, 0.0, b, 0.0,
//                          0.0, 0.0, 0.0, b;
//                     X = P * X * P.transpose();


//                     // 楕円Bを定義
//                     Matrix4f Delta, Q;
//                     Vector4f q_list[2];
//                     FTYPE qx = cos(theta1);
//                     FTYPE qy = sin(theta1) * cos(theta2);
//                     FTYPE qz = sin(theta1) * sin(theta2) * cos(theta3);
//                     FTYPE qw = sin(theta1) * sin(theta2) * sin(theta3);
//                     Q << qx,-qy,-qz,-qw,
//                          qy, qx, qw,-qz,
//                          qz,-qw, qx, qy,
//                          qw, qz,-qy, qx;
//                     FTYPE alpha = r * (delta*delta);
//                     FTYPE beta = r * sqrt(1 - (1-delta*delta)*(1-delta*delta));
//                     // std::cout << "α = " << alpha << std::endl;
//                     // std::cout << "β = " << beta << std::endl;
//                     Delta << alpha*alpha, 0.0, 0.0, 0.0,
//                              0.0, beta*beta, 0.0, 0.0,
//                              0.0, 0.0, beta*beta, 0.0,
//                              0.0, 0.0, 0.0, beta*beta;
//                     Delta = Q * Delta * Q.transpose();
//                     q_list[0] << r * (1 - delta*delta), 0, 0, 0;
//                     q_list[0] = Q * q_list[0];
//                     for(int j = 0; j < 4; j++) q_list[1](j) = q_list[0](j) + inv_sqrt2;

//                     // 基底行列の定義
//                     Matrix4f Y;
//                     Y << alpha, 0.0, 0.0, 0.0,
//                           0.0, beta, 0.0, 0.0,
//                           0.0, 0.0, beta, 0.0,
//                           0.0, 0.0, 0.0, beta;
//                     Y = Q * Y * Q.transpose();

//                     // ユニモジュラ行列の計算と作用
//                     start = std::chrono::system_clock::now();
//                     UnimodularMatrix<4> G_T = LLL_ZRoot2<4>(X, Y, (FTYPE)0.90);
//                     // G_T = UnimodularMatrix<4>::Identity();
//                     UnimodularMatrix<4> G = G_T.transpose();
//                     UnimodularMatrix<4> G_dot = conj(G);
//                     UnimodularMatrix<4> G_inv = inverse_UnimodularMatrix(G);
                    
//                     Matrix4f G_FTYPE = ZRoot2_to_FTYPE(G);
//                     Matrix4f G_dot_FTYPE = ZRoot2_to_FTYPE(G_dot);
//                     Matrix4f G_inv_FTYPE = ZRoot2_to_FTYPE(G_inv);

//                     for(int i = 0; i < 2; i++){
//                         p_list[i] = G_FTYPE * p_list[i];
//                         q_list[i] = G_dot_FTYPE * q_list[i];
//                     }
//                     D = G_FTYPE * D * G_FTYPE.transpose();
//                     Delta = G_dot_FTYPE * Delta * G_dot_FTYPE.transpose();

//                     end = std::chrono::system_clock::now();
//                     time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
//                     // std::cout << "time1 " << time << "[ms]" << std::endl;
                    
                    
//                     start = std::chrono::system_clock::now();
//                     for(int u_omega = 0; u_omega < 2; u_omega++){
//                         for(int t_omega = 0; t_omega < 2; t_omega++){
//                             Vector4f p = p_list[u_omega];
//                             Vector4f q = q_list[t_omega];
//                             std::array<std::pair<FTYPE, FTYPE>, 4> A, B;
//                             for(int j = 0; j < 4; j++){
//                                 A[j] = {p(j) - sqrt(D(j,j)), p(j) + sqrt(D(j,j))};
//                                 B[j] = {q(j) - sqrt(Delta(j,j)), q(j) + sqrt(Delta(j,j))};
//                                 // std::cout << A[j].second - A[j].first << std::endl;
//                                 // std::cout << B[j].second - B[j].first << std::endl;
//                                 if((A[j].second - A[j].first) * (B[j].second - B[j].first) > 1000){
//                                     std::cout << "delta " << (A[j].second - A[j].first) * (B[j].second - B[j].first) << std::endl;
//                                 }
//                                 // std::cout << std::endl;
//                             }


//                             std::vector<ZRoot2> X = one_dim_grid_problem(A[0].first, A[0].second, 
//                                                                         B[0].first, B[0].second);

//                             std::vector<ZRoot2> Y = one_dim_grid_problem(A[1].first, A[1].second, 
//                                                                         B[1].first, B[1].second);

//                             std::vector<ZRoot2> Z = one_dim_grid_problem(A[2].first, A[2].second, 
//                                                                         B[2].first, B[2].second);

//                             std::vector<ZRoot2> W = one_dim_grid_problem(A[3].first, A[3].second, 
//                                                                         B[3].first, B[3].second);


//                             // std::cout << X.size() << " " << Y.size() << " " << Z.size() << " " << W.size() << std::endl;
//                             // std::cout << "解の個数 : " << X.size() * Y.size() * Z.size() * W.size() << std::endl;
//                             if(X.size() * Y.size() * Z.size() * W.size() > 100000){
//                                 std::cout << "解の個数 : " << X.size() * Y.size() * Z.size() * W.size() << std::endl;
//                             }

//                             // i=0の場合は(√2x)^2 + (√2y)^2 + (√2z)^2 + (√2w)^2 = 2^{k+1}
//                             // i=1の場合は(√2x + 1)^2 + (√2y + 1)^2 + (√2z + 1)^2 + (√2w + 1)^2 = 2^{k+1}
//                             // を満たす(x,y,z,w)を探す。
//                             // for(auto &x : X) x *= ZRoot2(0, 1);
//                             // for(auto &y : Y) y *= ZRoot2(0, 1);
//                             // for(auto &z : Z) z *= ZRoot2(0, 1);
//                             // for(auto &w : W) w *= ZRoot2(0, 1);

                            
//                             ZRoot2 key = {((ITYPE)1<<(k)), 0};
//                             for(auto &x : X){
//                                 for(auto &y : Y){
//                                     for(auto &z : Z){
//                                         for(auto &w : W){
//                                             Eigen::Matrix<ZRoot2, 4, 1> xyzw;
//                                             xyzw << x, y, z, w;
//                                             xyzw = G_inv * xyzw;
//                                             ITYPE a0 = xyzw(0).a;
//                                             ITYPE b0 = xyzw(0).b;
//                                             ITYPE a1 = xyzw(1).a;
//                                             ITYPE b1 = xyzw(1).b;
//                                             ITYPE a2 = xyzw(2).a;
//                                             ITYPE b2 = xyzw(2).b;
//                                             ITYPE a3 = xyzw(3).a;
//                                             ITYPE b3 = xyzw(3).b;

//                                             ZOmega u,t;
//                                             if(u_omega == 0) u = {-b0+b1, a1, b0+b1, a0};
//                                             else u = {-b0+b1, a1, b0+b1+1, a0};
                                                
//                                             if(t_omega == 0) t = {-b2+b3, a3, b2+b3, a2};
//                                             else t = {-b2+b3, a3, b2+b3+1, a2};
                
//                                             if(u.norm() + t.norm() == key){
//                                                 U2_ZOmega V(u, t, l, k);
//                                                 if(distance(U, to_quaternion(V)) < eps){
// #pragma omp critical
//                                                     V_list.push_back(ExactSynthesis(V));
//                                                 }
//                                             }
//                                         }
//                                     }
//                                 }
//                             }
//                         }
//                         end = std::chrono::system_clock::now();
//                         time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
//                         // std::cout << "time2 " << time << "[ms]" << std::endl;
//                         }
//                 }
//             }
//         }
        
//         if(!V_list.empty()){
//             std::cout << "Start Exact Synthesis" << std::endl;
//             auto comp = [](const std::string& x, const std::string &y){
//                 return std::count(x.begin(), x.end(), 'T') < std::count(y.begin(), y.end(), 'T');
//             };
//             return *min_element(V_list.begin(), V_list.end(), comp);
//         }
//     }
//     return "Failure";
// }



// std::string su2compiler(quaternion U, FTYPE eps, int l){
//     Matrix4f P;
//     P << U.a,-U.b,-U.c,-U.d,
//          U.b, U.a, U.d,-U.c,
//          U.c,-U.d, U.a, U.b,
//          U.d, U.c,-U.b, U.a;
    
//     Matrix8f SIGMA, inv_SIGMA;
//     SIGMA <<    1.0, sqrt2,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
//                 0.0,   0.0,  1.0, sqrt2,   0.0,   0.0,   0.0,   0.0,
//                 0.0,   0.0,  0.0,   0.0,   1.0, sqrt2,   0.0,   0.0,
//                 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   1.0, sqrt2,
//                 1.0,-sqrt2,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
//                 0.0,   0.0,  1.0,-sqrt2,   0.0,   0.0,   0.0,   0.0,
//                 0.0,   0.0,  0.0,   0.0,   1.0,-sqrt2,   0.0,   0.0,
//                 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   1.0,-sqrt2;
//     inv_SIGMA = SIGMA.transpose();
//     inv_SIGMA.row(0) /= 2.0;
//     inv_SIGMA.row(2) /= 2.0;
//     inv_SIGMA.row(4) /= 2.0;
//     inv_SIGMA.row(6) /= 2.0;
//     inv_SIGMA.row(1) /= 4.0;
//     inv_SIGMA.row(3) /= 4.0;
//     inv_SIGMA.row(5) /= 4.0;
//     inv_SIGMA.row(7) /= 4.0;


//     for(int k = 0; k < 100; k++){
//         std::cout << "k = " << k << std::endl;
//         FTYPE r = pow(sqrt2, (FTYPE)k);
//         std::vector<std::string> V_list;

//         FTYPE delta = min((FTYPE)1.0, 3.0*pow(r, -(FTYPE)1.6) / eps);
//         std::cout << "delta = " << delta << std::endl;
//         std::cout << "a = " << sqrt(8.0/(PI*PI) * (1 - pow(1-delta*delta, 1.0/(FTYPE)3.0))) << std::endl;
//         int N = (int)ceil(1.0/ sqrt(8.0/(PI*PI) * (1 - pow(1-delta*delta, 1.0/(FTYPE)3.0))));
//         std::cout << "N = " << N << std::endl;

//         std::vector<FTYPE> theta1_list(N), theta2_list(N), theta3_list(2*N);
//         for(int i = 0; i < N; i++) theta1_list[i] = PI * ((FTYPE)i / N);
//         for(int i = 0; i < N; i++) theta2_list[i] = PI * ((FTYPE)i / N);
//         for(int i = 0; i < 2*N; i++) theta3_list[i] = PI * ((FTYPE)i / N);

// #pragma omp parallel for collapse(3)
//         for(auto &theta1 : theta1_list){
//             for(auto &theta2 : theta2_list){
//                 for(auto &theta3 : theta3_list){
//                     FTYPE qx = cos(theta1);
//                     FTYPE qy = sin(theta1) * cos(theta2);
//                     FTYPE qz = sin(theta1) * sin(theta2) * cos(theta3);
//                     FTYPE qw = sin(theta1) * sin(theta2) * sin(theta3);
//                     FTYPE alpha = r * (delta*delta);
//                     FTYPE beta = r * sqrt(1 - (1-delta*delta)*(1-delta*delta));

//                     start = std::chrono::system_clock::now();
//                     // std::cout << "start" << std::endl;

//                     // X基底を定義
//                     Matrix4f X;
//                     FTYPE a = r * (1 - sqrt(1 - eps*eps));
//                     FTYPE b = r * eps;
//                     X << a, 0.0, 0.0, 0.0,
//                          0.0, b, 0.0, 0.0,
//                          0.0, 0.0, b, 0.0,
//                          0.0, 0.0, 0.0, b;
//                     X = P * X * P.transpose();

//                     // Y基底を定義
//                     Matrix4f Y, Q;
//                     Y << alpha, 0.0, 0.0, 0.0,
//                           0.0, beta, 0.0, 0.0,
//                           0.0, 0.0, beta, 0.0,
//                           0.0, 0.0, 0.0, beta;
//                     Q << qx,-qy,-qz,-qw,
//                          qy, qx, qw,-qz,
//                          qz,-qw, qx, qy,
//                          qw, qz,-qy, qx;
//                     Y = Q * Y * Q.transpose();

//                     // B基底を定義
//                     Matrix8f B;
//                     for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) B(i,j) = X(i,j);
//                     for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) B(i+4,j+4) = Y(i,j);
//                     B = B * SIGMA;


//                     // ユニモジュラ行列を計算
//                     UnimodularMatrix<8> G_T = LLL<8>(B, 0.5);
//                     UnimodularMatrix<8> G = G_T.transpose();
//                     UnimodularMatrix<8> inv_G = inverse_UnimodularMatrix(G);
//                     Matrix8f G_FTYPE;
//                     for(int i = 0; i < 8; i++) for(int j = 0; j < 8; j++) G_FTYPE(i,j) = (FTYPE)G(i,j);

//                     // 楕円行列の計算
//                     // x^T D^{-1} x ≦ 1なので、G^Tになってることに注意
//                     Matrix8f D;
//                     D = G_FTYPE * B.transpose() * B * G_FTYPE.transpose();


//                     // 中心点の計算
//                     Vector4f p_list[2];
//                     p_list[0] << r * sqrt(1 - eps*eps), 0, 0, 0;
//                     p_list[0] << P * p_list[0];
//                     for(int j = 0; j < 4; j++) p_list[1](j) = p_list[0](j) - inv_sqrt2;

//                     Vector4f q_list[2];
//                     q_list[0] << r * (1 - delta*delta), 0, 0, 0;
//                     q_list[0] = Q * q_list[0];
//                     for(int j = 0; j < 4; j++) q_list[1](j) = q_list[0](j) + inv_sqrt2;

//                     Vector8f pq_list[2];
//                     for(int i = 0; i < 2; i++){
//                         for(int j = 0; j < 4; j++) pq_list[i](j) = p_list[i](j);
//                         for(int j = 0; j < 4; j++) pq_list[i](j+4) = q_list[i](j);
//                         pq_list[i] = inv_SIGMA * pq_list[i];
//                         pq_list[i] = G_FTYPE * pq_list[i];
//                     }
                
//                     end = std::chrono::system_clock::now();
//                     time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
//                     // std::cout << "time1 " << time << "[ms]" << std::endl;

//                     // 各中心点に対して格子点列挙
//                     ZRoot2 key = {(ITYPE)1<<(k), 0};
//                     for(int i = 0; i < 2; i++){
//                         Vector8f pq = pq_list[i];
//                         ITYPE MIN[8];
//                         ITYPE MAX[8];
//                         ITYPE SIZE = 1;
//                         for(int j = 0; j < 8; j++){
//                             MIN[j] = (ITYPE)ceil(pq(j) - sqrt(D(j,j)));
//                             MAX[j] = (ITYPE)floor(pq(j) + sqrt(D(j,j)));
//                             SIZE *= (MAX[i] - MIN[i]);
//                         }
//                         // std::cout << "解の個数 " << SIZE << std::endl;

//                         for(ITYPE a0 = MIN[0]; a0 <= MAX[0]; a0++){
//                         for(ITYPE b0 = MIN[1]; b0 <= MAX[1]; b0++){
//                         for(ITYPE a1 = MIN[2]; a1 <= MAX[2]; a1++){
//                         for(ITYPE b1 = MIN[3]; b1 <= MAX[3]; b1++){  
//                         for(ITYPE a2 = MIN[4]; a2 <= MAX[4]; a2++){
//                         for(ITYPE b2 = MIN[5]; b2 <= MAX[5]; b2++){
//                         for(ITYPE a3 = MIN[6]; a3 <= MAX[6]; a3++){
//                         for(ITYPE b3 = MIN[7]; b3 <= MAX[7]; b3++){
//                             Eigen::Vector<ITYPE, 8> x;
//                             x << a0, b0, a1, b1, a2, b2, a3, b3;
//                             x = inv_G * x;
//                             ZOmega u,t;
//                             if(i == 0){
//                                 u = {-b0+b1, a1, b0+b1, a0};
//                                 t = {-b2+b3, a3, b2+b3, a2};
//                             }else{
//                                 u = {-b0+b1, a1, b0+b1, a0+1};
//                                 t = {-b2+b3, a3, b2+b3, a2+1};
//                             }
//                             if(u.norm() + t.norm() == key){
//                                 U2_ZOmega V(u, t, l, k);
//                                 if(distance(U, to_quaternion(V)) < eps){
//     #pragma omp critical
//                                     V_list.push_back(ExactSynthesis(V));
//                                 }
//                             }
//                         }
//                         }
//                         }
//                         }
//                         }
//                         }
//                         }
//                         }
//                     }
//                     end = std::chrono::system_clock::now();
//                     time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
//                     // std::cout << "time2 " << time << "[ms]" << std::endl;
//                 }
//             }
//         }
        
//         if(!V_list.empty()){
//             std::cout << "Start Exact Synthesis" << std::endl;
//             auto comp = [](const std::string& x, const std::string &y){
//                 return std::count(x.begin(), x.end(), 'T') < std::count(y.begin(), y.end(), 'T');
//             };
//             return *min_element(V_list.begin(), V_list.end(), comp);
//         }
//     }
//     return "Failure";
// }




}