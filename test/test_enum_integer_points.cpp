#include <gtest/gtest.h>

#include <iostream>
#include <algorithm>
#include <random>
#include <Eigen/Dense>
#include <Eigen/QR>

#include "type.hpp"
#include "enum_integer_points.hpp"

using namespace std;
using namespace Eigen;
using namespace su2compiler;

MatrixXR randomOrthogonalMatrix(int n, unsigned int seed=1234)
{
    mt19937 gen(seed);
    uniform_real_distribution<> dist(-1.0, 1.0);

    // ランダム行列生成（成分は一様分布 [0,1)）
    MatrixXd A(n,n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A(i, j) = dist(gen);

    // QR分解
    HouseholderQR<MatrixXd> qr(A);
    MatrixXd Q = qr.householderQ();

    MatrixXR Q_Real(n,n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++) Q_Real(i,j) = Q(i,j);
    }
    return Q_Real;
}

void test_GSO(){
    int m = 8;
    int n = 8;
    MatrixXR B = MatrixXR::Random(m,n);

    MatrixXR B_orth = GSO(B);

    Real tol = 1e-10;
    MatrixXR D1 = B.transpose() * B;
    MatrixXR D2 = B_orth.transpose() * B_orth;
    assert(abs(D1.determinant() - D2.determinant()) < tol);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i != j) assert(abs(D2(i,j)) < tol);
        }
    }
}

void test_LLL()
{
    int n = 8;
    // 基底行列を定義
    MatrixXR P = randomOrthogonalMatrix(n);   // 固有ベクトル
    Real e1 = 1e-10;   // 固有値
    Real e2 = 1;
    MatrixXR B = MatrixXR::Zero(n,n);
    for(int i = 0; i < n; i++){
        if(i == 0) B(i,i) = e1;
        else       B(i,i) = e2;
    }
    B = B * P;

    auto [U, U_inv] = LLL(B, 0.75);
    MatrixXR U_Real = U.cast<Real>();
    MatrixXR B_reduce = B * U_Real;

    Real delta = 1.0;
    for(int i = 0; i < n; i++){
        delta *= B_reduce.col(i).norm();
    }
    delta /= B_reduce.determinant();

    assert(1.0 <= delta && delta <= pow((Real)2, (Real)(n*(n-1)) / 4.0));

    MatrixXI UU_inv = U * U_inv;
    assert((UU_inv - MatrixXI::Identity(n,n)).cwiseAbs().maxCoeff() == 0);
}

void test_cholesky()
{
    int n = 8;
    // 正定置行列を定義
    MatrixXR P = randomOrthogonalMatrix(n);   // 固有ベクトル
    Real e1 = 1e-10;   // 固有値
    Real e2 = 1;
    MatrixXR A = MatrixXR::Zero(n,n);  
    for(int i = 0; i < n; i++){
        if(i == 0) A(i,i) = e1;
        else       A(i,i) = e2;
    }
    A = P.transpose() * A * P;

    MatrixXR L = cholesky(A);
    MatrixXR LLT = L * L.transpose();

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            assert(abs(LLT(i,j) - A(i,j)) < Real(1e-10));
        }
    }
}

void test_EnumIntegerPoints()
{
    int n = 8;
    // 基底行列を定義
    MatrixXR P = randomOrthogonalMatrix(n);   // 固有ベクトル
    Real e1 = 1e-1;   // 固有値
    Real e2 = 1e1;
    MatrixXR Q_inv(n,n);
    for(int i = 0; i < n; i++){
        if(i < n/2) Q_inv(i,i) = e1;
        else        Q_inv(i,i) = e2;
    }
    Q_inv = P.transpose() * Q_inv * P;
    MatrixXR Q = Q_inv.inverse();

    VectorXR p = VectorXR::Random(n);
    Real c = 1.0;
    
    std::vector<VectorXI> X = EnumIntegerPoints(Q, p, c);

    // Brute-force的に整数点探索
    std::vector<VectorXI> Y;
    VectorXI xvec(n);
    std::function<void(int)> BF_EnumIntegerPoints;
    BF_EnumIntegerPoints = [&Q_inv, &Q, &p, &c, &xvec, &Y, &BF_EnumIntegerPoints](int i)
    {
        if(i == Q.rows()){
            VectorXR xvec_Real = xvec.cast<Real>(); 
            if((xvec_Real - p).dot(Q * (xvec_Real - p)) <= c) Y.push_back(xvec);
        }else{
            Integer xi_min = ceil_mpreal(p(i) - sqrt(Q_inv(i,i) * c));
            Integer xi_max = floor_mpreal(p(i) + sqrt(Q_inv(i,i) * c));
            for(Integer xi = xi_min; xi <= xi_max; xi++){
                xvec(i) = xi;
                BF_EnumIntegerPoints(i+1);
            }
        }
    };
    BF_EnumIntegerPoints(0);

    auto comp = [](VectorXI x, VectorXI y){
        for(int i = 0; i < x.rows(); i++){
            if(x(i) < y(i)) return true;
            if(x(i) > y(i)) return false;
        }
        return false;
    };
    std::sort(X.begin(), X.end(), comp);
    std::sort(Y.begin(), Y.end(), comp);

    std::cout << "|X|=" << X.size() << std::endl;
    for(auto x : X) cout << x.transpose() << endl;
    std::cout << "|Y|=" << Y.size() << std::endl;
    for(auto y : Y) cout << y.transpose() << endl;

    assert(X == Y);
}


// int main()
// {
//     // test_GSO();
//     // test_LLL();
//     // test_cholesky();
//     test_EnumIntegerPoints();
// }