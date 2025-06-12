#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <Eigen/Dense>
#include <Eigen/QR>

#include "core/type.hpp"
#include "math/functions.hpp"
#include "math/lattice.hpp"

using namespace std;
using namespace su2compiler;


MatrixR RandomOrthogonalMatrix(int N, unsigned int seed=1234)
{
    mt19937 gen(seed);
    uniform_real_distribution<> dist(-1.0, 1.0);

    // ランダム行列生成（成分は一様分布 [0,1)）
    Eigen::MatrixXd A(N,N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            A(i, j) = dist(gen);

    // QR分解
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
    Eigen::MatrixXd Q = qr.householderQ();

    MatrixR Q_Real = Q.cast<Real>();
    return Q_Real;
}


TEST(LatticeTest, LLL)
{
    Real::set_default_prec(256);

    int N = 8;
    // 基底行列を定義
    MatrixR P = RandomOrthogonalMatrix(N);   // 固有ベクトル
    Real e1 = 1e-10;   // 固有値
    Real e2 = 1;
    MatrixR B = MatrixR::Zero(N,N);
    for(int i = 0; i < N; i++){
        if(i == 0) B(i,i) = e1;
        else       B(i,i) = e2;
    }
    B = B * P;

    auto [U, U_inv] = math::lattice::LLL(B, 0.75);
    MatrixR U_Real = U.unaryExpr([](Integer z) {
        return type::Int_to_Real(z);
    });
    MatrixR B_reduce = B * U_Real;

    Real delta = 1.0;
    for(int i = 0; i < N; i++){
        delta *= B_reduce.col(i).norm();
    }
    delta /= B_reduce.determinant();

    EXPECT_TRUE(1.0 <= delta && delta <= pow(Real(2), (Real)(N*(N-1)) / 4.0));

    MatrixI UU_inv = U * U_inv;
    EXPECT_TRUE((UU_inv - MatrixI::Identity(N,N)).cwiseAbs().maxCoeff() == 0);
}


TEST(LatticeTest, EnumIntegerPoints)
{
    int N = 8;
    // 基底行列を定義
    MatrixR P = RandomOrthogonalMatrix(N);   // 固有ベクトル
    Real e1 = 1e-1;   // 固有値
    Real e2 = 1e1;
    MatrixR Q_inv = MatrixR::Zero(N,N);
    for(int i = 0; i < N; i++){
        if(i < N/2) Q_inv(i,i) = e1;
        else        Q_inv(i,i) = e2;
    }
    Q_inv = P.transpose() * Q_inv * P;
    MatrixR Q = Q_inv.inverse();

    VectorR p = VectorR::Random(N);
    Real c = 1.0;

    std::vector<VectorI> X = math::lattice::EnumIntegerPoints(Q, p, c);

    // Brute-force的に整数点探索
    std::vector<VectorI> Y;
    VectorI xvec(N);
    std::function<void(int)> BF_EnumIntegerPoints;
    BF_EnumIntegerPoints = [&](int i)
    {
        if(i == Q.rows()){
            VectorR xvec_Real = xvec.unaryExpr([](Integer z) {
                return type::Int_to_Real(z);
            }); 
            if((xvec_Real - p).dot(Q * (xvec_Real - p)) <= c) Y.push_back(xvec);
        }else{
            Integer xi_min = math::ceil(p(i) - sqrt(Q_inv(i,i) * c));
            Integer xi_max = math::floor(p(i) + sqrt(Q_inv(i,i) * c));
            for(Integer xi = xi_min; xi <= xi_max; xi++){
                xvec(i) = xi;
                BF_EnumIntegerPoints(i+1);
            }
        }
    };
    BF_EnumIntegerPoints(0);

    auto comp = [](VectorI x, VectorI y){
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

