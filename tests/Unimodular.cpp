#include <bits/stdc++.h>
#include <Eigen/Core>
#include "src/type.hpp"
#include "src/rings.hpp"
#include "src/quaternion.hpp"
#include "src/Unimodular.hpp"
#include <chrono>

using namespace std;
using namespace SU2_Compiler;


int main(){
    // set_random_unitary_seed(1234);
    quaternion U = random_unitary();
    FTYPE eps = 1e-6;
    FTYPE a = eps*eps;
    FTYPE b = sqrt(1.0 - (1.0 - a) * (1.0 - a));
    Matrix4f D, Delta, P;
    D << 1 / (a*a), 0.0, 0.0, 0.0,
         0.0, 1 / (b*b), 0.0, 0.0,
         0.0, 0.0, 1 / (b*b), 0.0,
         0.0, 0.0, 0.0, 1 / (b*b);
    P << U.a,-U.b,-U.c,-U.d,
         U.b, U.a, U.d,-U.c,
         U.c,-U.d, U.a, U.b,
         U.d, U.c,-U.b, U.a;
    D = P * D * P.transpose();
    // D /= sqrt(sqrt(1 / (a*a * b*b * b*b * b*b)));
    Delta = Matrix4f::Identity();

    Matrix4f X,Y;
    X << 1 / a, 0.0, 0.0, 0.0,
         0.0, 1 / b, 0.0, 0.0,
         0.0, 0.0, 1 / b, 0.0,
         0.0, 0.0, 0.0, 1 / b;
    X = P * X * P.transpose();
    Y = Matrix4f::Identity();
    
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    UnimodularMatrix<4> G = LLL_ZRoot2<4>(X, Y, (FTYPE)0.75);
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    cout << "time : " << time << "[ms]" << endl; 
    UnimodularMatrix<4> G_dot = conj(G);

    cout << "G = \n";
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++) cout << G(i,j) << " ";
        cout << endl;
    }

    Matrix4f G_FTYPE = ZRoot2_to_FTYPE<4>(G);
    Matrix4f G_dot_FTYPE = ZRoot2_to_FTYPE<4>(conj<4>(G));

    cout << endl;
    cout << "D = \n" << D << endl;
    cout << "Δ = \n" << Delta << endl;
    cout << "D'= \n" << G_FTYPE.transpose() * D * G_FTYPE << endl;
    cout << "Δ'= \n" << G_dot_FTYPE.transpose() * Delta * G_dot_FTYPE << endl;
    cout << determinant(G) << endl;
    cout << determinant(G_dot) << endl; 
    cout << determinant(D) << endl;
    
    UnimodularMatrix<4> G_inv = inverse_UnimodularMatrix(G);
    Matrix4f G_inv_FTYPE = ZRoot2_to_FTYPE<4>(G_inv);
    cout << G_FTYPE * G_inv_FTYPE << endl;
    UnimodularMatrix<4> tmp;
    tmp = G * G_inv;
    cout << tmp << endl;
}