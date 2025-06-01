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


template <typename Real>
class LatticeTest : public ::testing::Test {
protected:
    void SetUp() override {
        if constexpr(std::is_same_v<Real, mpfr::mpreal>) {
            mpfr::mpreal::set_default_prec(256);
        }
    }
};
using RealTypes = ::testing::Types<REAL_SCALAR_TYPE_LIST>;
TYPED_TEST_SUITE(LatticeTest, RealTypes);


template <typename Real>
Eigen::MatrixX<Real> RandomOrthogonalMatrix(int n, unsigned int seed=1234)
{
    mt19937 gen(seed);
    uniform_real_distribution<> dist(-1.0, 1.0);

    // ランダム行列生成（成分は一様分布 [0,1)）
    Eigen::MatrixXd A(n,n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A(i, j) = dist(gen);

    // QR分解
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
    Eigen::MatrixXd Q = qr.householderQ();

    Eigen::MatrixX<Real> Q_Real = Q.cast<Real>();
    return Q_Real;
}


TYPED_TEST(LatticeTest, LLL)
{
    using Real = TypeParam;

    int n = 8;

    // 基底行列を定義
    Eigen::MatrixX<Real> P = RandomOrthogonalMatrix<Real>(n);   // 固有ベクトル
    Real e1 = 1e-10;   // 固有値
    Real e2 = 1;
    Eigen::MatrixX<Real> B = Eigen::MatrixX<Real>::Zero(n,n);
    for(int i = 0; i < n; i++){
        if(i == 0) B(i,i) = e1;
        else       B(i,i) = e2;
    }
    B = B * P;

    auto [U, U_inv] = math::lattice::LLL<Real>(B, 0.75);
    Eigen::MatrixX<Real> U_Real = U.unaryExpr([](Integer z) {
        return type::Int_to_Real<Real>(z);
    });
    Eigen::MatrixX<Real> B_reduce = B * U_Real;

    Real delta = 1.0;
    for(int i = 0; i < n; i++){
        delta *= B_reduce.col(i).norm();
    }
    delta /= B_reduce.determinant();

    EXPECT_TRUE(1.0 <= delta && delta <= pow((Real)2, (Real)(n*(n-1)) / 4.0));

    Eigen::MatrixX<Integer> UU_inv = U * U_inv;
    EXPECT_TRUE((UU_inv - Eigen::MatrixX<Integer>::Identity(n,n)).cwiseAbs().maxCoeff() == 0);
}


TYPED_TEST(LatticeTest, EnumIntegerPoints)
{
    using Real = TypeParam;

    int n = 8;
    // 基底行列を定義
    Eigen::MatrixX<Real> P = RandomOrthogonalMatrix<Real>(n);   // 固有ベクトル
    Real e1 = 1e-1;   // 固有値
    Real e2 = 1e1;
    Eigen::MatrixX<Real> Q_inv = Eigen::MatrixX<Real>::Zero(n,n);
    for(int i = 0; i < n; i++){
        if(i < n/2) Q_inv(i,i) = e1;
        else        Q_inv(i,i) = e2;
    }
    Q_inv = P.transpose() * Q_inv * P;
    Eigen::MatrixX<Real> Q = Q_inv.inverse();

    Eigen::VectorX<Real> p = Eigen::VectorX<Real>::Random(n);
    Real c = 1.0;

    // std::vector<Eigen::VectorX<Integer>> X;
    std::vector<Eigen::VectorX<Integer>> X = math::lattice::EnumIntegerPoints(Q, p, c);

    // Brute-force的に整数点探索
    std::vector<Eigen::VectorX<Integer>> Y;
    Eigen::VectorX<Integer> xvec(n);
    std::function<void(int)> BF_EnumIntegerPoints;
    BF_EnumIntegerPoints = [&](int i)
    {
        if(i == Q.rows()){
            Eigen::VectorX<Real> xvec_Real = xvec.unaryExpr([](Integer z) {
                return type::Int_to_Real<Real>(z);
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

    auto comp = [](Eigen::VectorX<Integer> x, Eigen::VectorX<Integer> y){
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

