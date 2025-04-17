#include <bits/stdc++.h>
#include <Eigen/Core>
#include "src/type.hpp"
#include "src/quaternion.hpp"
#include "src/Unimodular_Z.hpp"
#include <chrono>

using namespace std;
using namespace SU2_Compiler;

template <int M, int N>
FTYPE cost(Eigen::Matrix<FTYPE, M, N> B){
    Eigen::Matrix<FTYPE, N, N> D = B.transpose() * B;
    FTYPE det = determinant(D);
    FTYPE Prod = 1.0;
    for(int i = 0; i < N; i++) Prod *= D(i,i);
    return sqrt(Prod / det);
}


int main(){
    set_random_unitary_seed(1234);
    quaternion U = random_unitary();
    FTYPE eps = 1e-5;
    FTYPE r = pow(sqrt2, 25);
     
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
    // P << U.a,-U.b,-U.c,-U.d, 0, 0, 0, 0,
    //      U.b, U.a, U.d,-U.c, 0, 0, 0, 0,
    //      U.c,-U.d, U.a, U.b, 0, 0, 0, 0,
    //      U.d, U.c,-U.b, U.a, 0, 0, 0, 0,
    //      0, 0, 0, 0, U.a,-U.b,-U.c,-U.d,
    //      0, 0, 0, 0, U.b, U.a, U.d,-U.c, 
    //      0, 0, 0, 0, U.c,-U.d, U.a, U.b,
    //      0, 0, 0, 0, U.d, U.c,-U.b, U.a;

    Matrix8f Sigma, Sigma_inv;
    // Sigma << 1, 1/sqrt2, 0,-1/sqrt2, 0, 0, 0, 0,
    //          0, 1/sqrt2, 1, 1/sqrt2, 0, 0, 0, 0,
    //          0, 0, 0, 0, 1, 1/sqrt2, 0,-1/sqrt2,
    //          0, 0, 0, 0, 0, 1/sqrt2, 1, 1/sqrt2,
    //          1,-1/sqrt2, 0, 1/sqrt2, 0, 0, 0, 0,
    //          0,-1/sqrt2, 1,-1/sqrt2, 0, 0, 0, 0,
    //          0, 0, 0, 0, 1,-1/sqrt2, 0, 1/sqrt2,
    //          0, 0, 0, 0, 0,-1/sqrt2, 1,-1/sqrt2;
    Sigma << 1, 0, 0, 0, 1/sqrt2,-1/sqrt2, 0, 0,
             0, 1, 0, 0, 1/sqrt2, 1/sqrt2, 0, 0,
             0, 0, 1, 0, 0, 0, 1/sqrt2,-1/sqrt2,
             0, 0, 0, 1, 0, 0, 1/sqrt2, 1/sqrt2,
             1, 0, 0, 0,-1/sqrt2, 1/sqrt2, 0, 0, 
             0, 1, 0, 0,-1/sqrt2,-1/sqrt2, 0, 0,
             0, 0, 1, 0, 0, 0,-1/sqrt2, 1/sqrt2,
             0, 0, 0, 1, 0, 0,-1/sqrt2,-1/sqrt2;

    Sigma_inv = Sigma.transpose();
    Sigma_inv /= 2.0;

    Matrix8f H;
    H << 1, 0, 0, 0, 1, 0, 0, 0,
         0, 1, 0, 0, 0, 1, 0, 0,
         0, 0, 1, 0, 0, 0, 1, 0,
         0, 0, 0, 1, 0, 0, 0, 1,
         1, 0, 0, 0,-1, 0, 0, 0,
         0, 1, 0, 0, 0,-1, 0, 0,
         0, 0, 1, 0, 0, 0,-1, 0,
         0, 0, 0, 1, 0, 0, 0,-1;
    H /= sqrt2;


    Matrix8f B;
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
    // B = P * B * P.transpose();
    B = B * P.transpose();
    // B = B * H;
    B = B * Sigma_inv.transpose();

    Eigen::Matrix<FTYPE, 8, 4> B1, B2;
    for(int i = 0; i < 4; i++){
        B1.col(i) = B.col(i);
        B2.col(i) = B.col(i+4);
    }
    cout << "B = \n" << B << endl;
    
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    UnimodularMatrix<4> G1 = LLL<8,4>(B1, (FTYPE)0.5);
    UnimodularMatrix<4> G2 = LLL<8,4>(B2, (FTYPE)0.5);
    UnimodularMatrix<8> G = UnimodularMatrix<8>::Zero();
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            G(i,j) = G1(i,j);
            G(i+4,j+4) = G2(i,j);
        }
    }
    G = LLL<8,8>(B, 0.9);
    Eigen::Matrix<FTYPE, 4, 4> G1_FTYPE, G2_FTYPE;
    for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) G1_FTYPE(i,j) = (FTYPE)G1(i,j);
    for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) G2_FTYPE(i,j) = (FTYPE)G2(i,j);
    
    Eigen::Matrix<FTYPE, 8, 8> G_FTYPE;
    for(int i = 0; i < 8; i++) for(int j = 0; j < 8; j++) G_FTYPE(i,j) = (FTYPE)G(i,j);
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    cout << "time : " << time << "[ms]" << endl;

    cout << "G = \n";
    for(int i = 0; i < 8; i++){
        for(int j = 0; j < 8; j++) cout << G(i,j) << " ";
        cout << endl;
    }

    Matrix8f B_prime = B * G_FTYPE;
    cout << "B' = \n" << B_prime << endl;
    Eigen::Matrix<FTYPE, 8, 8> D, D_prime;
    D = B.transpose() * B;
    D_prime = B_prime.transpose() * B_prime;

    cout << endl;
    cout << "D = \n" << D << endl;
    cout << "D'= \n" << D_prime << endl;

    FTYPE det = determinant(D);
    FTYPE Prod = 1.0;
    for(int i = 0; i < 8; i++) Prod *= D(i,i);
    FTYPE det_prime = determinant(D_prime);
    FTYPE Prod_prime = 1.0;
    for(int i = 0; i < 8; i++) Prod_prime *= D_prime(i,i);

    cout << "det(D) = " << determinant(D) << endl;
    cout << "cost(D) = " << cost<8,8>(D) << endl;
    cout << "cost(D') = " << cost<8,8>(D_prime) << endl;
    cout << sqrt(Prod) << " " << sqrt(Prod_prime) << endl;

    Eigen::Matrix<FTYPE, 4, 4> D1, D2;
    D1 = B1.transpose() * B1;
    D2 = B2.transpose() * B2;
    cout << determinant(D) << endl;
    cout << determinant(D1) << endl;
    cout << determinant(D2) << endl;
    Eigen::Matrix<FTYPE, 8, 4> B1_prime, B2_prime;
    B1_prime = B1 * G1_FTYPE;
    B2_prime = B2 * G2_FTYPE;
    cout << cost<8,4>(B1_prime) << endl;
    cout << cost<8,4>(B2_prime) << endl;

    // for(int i = 0; i < 4; i++){
    //     for(int j = 0; j < 4; j++){
    //         D1(i,j) = D(i,j);
    //         D2(i,j) = D(i+4,j+4);
    //     }
    // }

    // cout << D << endl;
    // cout << determinant(D) << endl;
    // cout << endl;
    // cout << D1 << endl;
    // cout << determinant(D1) << endl;
    // cout << endl;
    // cout << D2 << endl;
    // cout << determinant(D2) << endl;
}