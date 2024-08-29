#include <bits/stdc++.h>
#include <Eigen/Core>
#include "src/type.hpp"
#include "src/rings.hpp"
#include "src/quaternion.hpp"
#include "src/grid_operator.hpp"

using namespace std;
using namespace SU2_Compiler;

// template <typename T>
// vector<vector<T>> getSubMatrix(vector<vector<T>>& A, int p, int q, int N){
//     vector<vector<T>> A_sub(N-1, vector<T>(N-1));
//     int i = 0, j = 0;
//     for(int row = 0; row < N; row++){
//         for(int col = 0; col < N; col++){
//             if(row == p || col == q) continue;

//             A_sub[i][j] = A[row][col];
//             j++;
//             if(j == N-1){
//                 j = 0;
//                 i++;
//             }
//         }
//     }
//     return A_sub;
// }

// template <typename T>
// T determinant(vector<vector<T>>& A, int N){
//     if(N == 1) return A[0][0];
    
//     T ret = 0;
//     for(int i = 0; i < N; i++){
//         vector<vector<T>> A_sub = getSubMatrix(A, 0, i, N);
//         if(i % 2) ret -= A[0][i] * determinant(A_sub, N-1);
//         else ret += A[0][i] * determinant(A_sub, N-1);
//     }
//     return ret;
// } 

// template <typename T, int N>
// T determinant(Eigen::Matrix<T, N, N>& A){
//     vector<vector<T>> A_vec(N, vector<T>(N));
//     for(int i = 0; i < N; i++){
//         for(int j = 0; j < N; j++) A_vec[i][j] = A(i,j);
//     }

//     return determinant(A_vec, N);
// }


int main(){
    // set_random_unitary_seed(1234);
    quaternion U = random_unitary();
    FTYPE eps = 1e-8;
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
    D /= sqrt(sqrt(1 / (a*a * b*b * b*b * b*b)));
    Delta = Matrix4f::Identity();

    Matrix4f X,Y;
    X << 1 / a, 0.0, 0.0, 0.0,
         0.0, 1 / b, 0.0, 0.0,
         0.0, 0.0, 1 / b, 0.0,
         0.0, 0.0, 0.0, 1 / b;
    X = P * X * P.transpose();
    Y = Matrix4f::Identity();
    
    GridOp<4> G = LLL_ZRoot2<4>(X, Y, (FTYPE)0.99);
    GridOp<4> G_dot = conj(G);

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
    
    GridOp<4> G_inv = inverse_Special_GridOp(G);
    Matrix4f G_inv_FTYPE = ZRoot2_to_FTYPE<4>(G_inv);
    cout << G_FTYPE * G_inv_FTYPE << endl;
    GridOp<4> tmp;
    tmp = G * G_inv;
    cout << tmp << endl;
}