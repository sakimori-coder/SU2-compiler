#pragma once

#include <string>
#include <Eigen/Core>

#include "type.hpp"
#include "ring/all.hpp"
#include "su2.hpp"

namespace su2compiler {
namespace clifford_t {

using Matrix3Zroot2 = Eigen::Matrix<ring::Zroot2,3,3>;
using Matrix2Zzeta8 = Eigen::Matrix<ring::Zzeta8,2,2>;

struct SO3Droot2 {
    Matrix3Zroot2 mat;
    int k;

    SO3Droot2() noexcept {
        mat << ring::Zroot2(1), ring::Zroot2(0), ring::Zroot2(0),
               ring::Zroot2(0), ring::Zroot2(1), ring::Zroot2(0),
               ring::Zroot2(0), ring::Zroot2(0), ring::Zroot2(1);
    }
    SO3Droot2(Matrix3Zroot2 _mat, int _k);

    [[nodiscard]] SO3Droot2 transpose() const { return {mat.transpose(), k}; }

    void reduce();
    SO3Droot2& operator*=(const SO3Droot2& r);
};

[[nodiscard]] inline SO3Droot2 operator*(SO3Droot2 lhs,const SO3Droot2& rhs) noexcept { return lhs *= rhs; }




struct U2Dzeta8 {
    Matrix2Zzeta8 mat;
    UINT k;

    U2Dzeta8() noexcept {
        mat << ring::Zzeta8(1), ring::Zzeta8(0),
               ring::Zzeta8(0), ring::Zzeta8(1);
        k = 0;
    }
    U2Dzeta8(Matrix2Zzeta8 _mat, UINT _k);
    U2Dzeta8(ring::Zzeta8j q, int l, UINT _k) : U2Dzeta8(q.u, q.t, l, _k) {}
    U2Dzeta8(ring::Zzeta8 u, ring::Zzeta8 t, int l, UINT _k);
    U2Dzeta8(std::string sequence);

    [[nodiscard]] U2Dzeta8 adjoint() const {
        Matrix2Zzeta8 mat_adj;
        mat_adj << mat(0,0).conj_complex(), mat(1,0).conj_complex(),
                   mat(0,1).conj_complex(), mat(1,1).conj_complex();
        return {mat_adj, k};
    }

    [[nodiscard]] Matrix2C to_Matrix2C() const {
        return (Matrix2C() << 
                mat(0,0).to_Complex(), mat(0,1).to_Complex(),
                mat(1,0).to_Complex(), mat(1,1).to_Complex()
               ).finished() / pow_ui(SQRT2(), k);
    }

    [[nodiscard]] SU2 to_SU2() const;

    [[nodiscard]] SO3Droot2 to_SO3Droot2() const;

    void reduce();
    U2Dzeta8& operator*=(const U2Dzeta8& r);
};

[[nodiscard]] inline U2Dzeta8 operator*(U2Dzeta8 lhs,const U2Dzeta8& rhs) noexcept { return lhs *= rhs; }



std::string exact_synthesis(U2Dzeta8 U);
std::string exact_synthesis(SO3Droot2 U);

int get_Tcount(const U2Dzeta8& U);

}
}