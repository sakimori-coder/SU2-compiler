#include <bits/stdc++.h>
#include <Eigen/Core>

#include "src/CPU/linalg.hpp"

using namespace std;

int main(){
    srand((unsigned int) time(0));

    const int N = 20;
    const int M = 40;
    Eigen::MatrixXd A(N, M);
    Eigen::VectorXd b(N);
    A = Eigen::MatrixXd::Random(N, M);
    b = Eigen::VectorXd::Random(N);

    // cout << A << endl;
    // cout << b << endl;

    auto [V, inter] = solve_linear_system<double, N, M>(A, b);

    // cout << V << endl;
    // cout << inter << endl;
    
    Eigen::VectorXd x(M);
    Eigen::VectorXd coef(M-N);
    coef = Eigen::VectorXd::Random(M-N);   // ランダム係数
    x = inter;
    for(int i = 0; i < M-N; i++) x += coef(i) * V.col(i);

    Eigen::VectorXd diff(N);
    diff = A*x - b;
    cout << diff.norm() << endl;   // ||Ax-b||を確認

    // cout << V.transpose() * V << endl;   // 正規直交性を確認
}