#pragma once

#include <vector>
#include <Eigen/Core>

#include "core/type.hpp"
#include "core/su2.hpp"

namespace su2compiler
{

template <typename Real>
class MixSU2
{
public:
    Eigen::VectorX<Real> prob;
    std::vector<SU2<Real>> U;

    MixSU2(const Eigen::VectorX<Real>& _prob, 
            const std::vector<SU2<Real>>& _U) noexcept : prob(_prob), U(_U) {}
    MixSU2(const std::vector<SU2<Real>>& _U) noexcept : U(_U) {}

    void compute_optimal_prob(
            const SU2<Real>& targetV,
            Real lambda,
            Real betaStar,
            Real betaBar,
            Real gamma,
            Real epsilon1,
            Real epsilon2,
            int  MaxIteration,
            bool PUTPUT_HISTORY = false);

};


namespace sdp_tool
{

template <typename Real>
struct BlockMatrix
{
    Eigen::Matrix4<Real> blk1;
    Eigen::Matrix4<Real> blk2;
    Eigen::VectorX<Real> blk3;
    Real                 blk4;

    inline void setIdentity(int N) {
        blk1 = Eigen::Matrix4<Real>::Identity();
        blk2 = Eigen::Matrix4<Real>::Identity();
        blk3 = Eigen::VectorX<Real>::Ones(N);
        blk4 = Real(1);
    }

    inline BlockMatrix<Real>& operator+=(const BlockMatrix<Real>& r) noexcept {
        blk1 += r.blk1;
        blk2 += r.blk2;
        blk3 += r.blk3;
        blk4 += r.blk4;
        return *this;
    }
    inline BlockMatrix<Real> operator+(const BlockMatrix<Real>& r) noexcept {
        return *this += r;
    }

    inline BlockMatrix<Real>& operator-=(const BlockMatrix<Real>& r) noexcept {
        blk1 -= r.blk1;
        blk2 -= r.blk2;
        blk3 -= r.blk3;
        blk4 -= r.blk4;
        return *this;
    }
    inline BlockMatrix<Real> operator-(const BlockMatrix<Real>& r) noexcept {
        return *this -= r;
    }

    inline BlockMatrix<Real>& operator*=(const BlockMatrix<Real>& r) noexcept {
        blk1 *= r.blk1;
        blk2 *= r.blk2;
        blk3 = blk3.cwiseProduct(r.blk3);
        blk4 *= blk4;
        return *this;
    }
    inline BlockMatrix<Real> operator*(Real r) noexcept {
        return *this *= r;
    }

    inline BlockMatrix<Real>& operator*=(const Real& r) noexcept {
        blk1 *= r;
        blk2 *= r;
        blk3 *= r;
        blk4 *= r;
        return *this;
    }
    inline BlockMatrix<Real> operator*(const BlockMatrix<Real>& r) noexcept {
        return *this *= r;
    }
};

template <typename Real>
Real HSinner(
    const BlockMatrix<Real>& A,
    const BlockMatrix<Real>& B) noexcept 
{
    Real ret(0);
    ret += A.blk1.cwiseProduct(B.blk1).sum();
    ret += A.blk2.cwiseProduct(B.blk2).sum();
    ret += A.blk3.cwiseProduct(B.blk3).sum();
    ret += A.blk4.cwiseProduct(B.blk4).sum();
    return ret;
}


}


}