#pragma once

#include <string>
#include <Eigen/Core>

#include "core/type.hpp"
#include "core/su2.hpp"
#include "ring/all.hpp"
#include "math/constants.hpp"
#include "math/functions.hpp"


namespace su2compiler::clifford_t::exact {

using Matrix3Zroot2 = Eigen::Matrix<ring::Zroot2,3,3>;
using Matrix2Zzeta8 = Eigen::Matrix<ring::Zzeta8,2,2>;

struct SO3Droot2 {
    Matrix3Zroot2 mat;
    Natural k;

    SO3Droot2() noexcept {
        mat << ring::Zroot2(1), ring::Zroot2(0), ring::Zroot2(0),
               ring::Zroot2(0), ring::Zroot2(1), ring::Zroot2(0),
               ring::Zroot2(0), ring::Zroot2(0), ring::Zroot2(1);
    }
    SO3Droot2(Matrix3Zroot2 _mat, Natural _k);

    [[nodiscard]] SO3Droot2 transpose() const { return {mat.transpose(), k}; }

    void reduce();
    SO3Droot2& operator*=(const SO3Droot2& r);
};

[[nodiscard]] inline SO3Droot2 operator*(SO3Droot2 lhs,const SO3Droot2& rhs) noexcept { return lhs *= rhs; }




struct U2Dzeta8 {
    Matrix2Zzeta8 mat;
    Natural k;

    U2Dzeta8() noexcept {
        mat << ring::Zzeta8(1), ring::Zzeta8(0),
               ring::Zzeta8(0), ring::Zzeta8(1);
        k = 0;
    }
    U2Dzeta8(Matrix2Zzeta8 _mat, Natural _k);
    U2Dzeta8(ring::Zzeta8j q, int l, Natural _k) : U2Dzeta8(q.u, q.t, l, _k) {}
    U2Dzeta8(ring::Zzeta8 u, ring::Zzeta8 t, int l, Natural _k);
    U2Dzeta8(std::string sequence);

    [[nodiscard]] U2Dzeta8 adjoint() const {
        Matrix2Zzeta8 mat_adj;
        mat_adj << mat(0,0).conj_complex(), mat(1,0).conj_complex(),
                   mat(0,1).conj_complex(), mat(1,1).conj_complex();
        return {mat_adj, k};
    }

    template <typename RealType>
    [[nodiscard]] Eigen::Matrix2<std::complex<RealType>> toMatrix2C() const {
        return (Eigen::Matrix2<std::complex<RealType>>() << 
                mat(0,0).toComplex<RealType>(), mat(0,1).toComplex<RealType>(),
                mat(1,0).toComplex<RealType>(), mat(1,1).toComplex<RealType>()
               ).finished() / math::pow_ui(math::SQRT2<RealType>(), k);
    }

    template <typename RealType>
    [[nodiscard]] SU2<RealType> toSU2() const;

    [[nodiscard]] SO3Droot2 toSO3Droot2() const;

    void reduce();
    U2Dzeta8& operator*=(const U2Dzeta8& r);
};

[[nodiscard]] inline U2Dzeta8 operator*(U2Dzeta8 lhs,const U2Dzeta8& rhs) noexcept { return lhs *= rhs; }



std::string synth(U2Dzeta8 U);
std::string synth(SO3Droot2 U);

Natural get_Tcount(const U2Dzeta8& U);

}