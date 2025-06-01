#include "math/lattice.hpp"

#include <utility>
#include <algorithm>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "core/type.hpp"
#include "math/functions.hpp"
#include "math/linalg.hpp"

#include "util/time_profiler.hpp"


namespace su2compiler::math::lattice{

template <typename Real>
std::pair<Eigen::MatrixX<Integer>, Eigen::MatrixX<Integer>> LLL(
        Eigen::MatrixX<Real> B,
        Real delta)

{
    auto& prof = util::Profiler::instance();
    prof.start("LLL");

    int M = B.rows();
    int N = B.cols();
    Eigen::MatrixX<Real> B_orth(M, N);
    Eigen::MatrixX<Real> mu(N, N);
    Eigen::VectorX<Real> norm2_B_orth(N);

    B_orth = math::linalg::GSO(B);
    for(int i = 0; i < N; i++) norm2_B_orth[i] = B_orth.col(i).squaredNorm();
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            mu(i, j) = B.col(i).dot(B_orth.col(j)) / B_orth.col(j).squaredNorm();
        }
    }

    Eigen::MatrixX<Integer> U = Eigen::MatrixX<Integer>::Identity(N, N);
    Eigen::MatrixX<Integer> U_inv = Eigen::MatrixX<Integer>::Identity(N, N);

    int k = 1;
    while(k < N){
        // std::cout << k << std::endl;

        for(int j = k-1; j >= 0; j--){
            if(math::abs(mu(k,j)) > 0.5) {
                Integer q = math::round(mu(k,j));
                Real q_real = type::Int_to_Real<Real>(q);
                B.col(k) -= q_real * B.col(j);
                U.col(k) -= q * U.col(j);
                U_inv.row(j) += q * U_inv.row(k);
                for(int l = 0; l <= j; l++) mu(k,l) -= q_real * mu(j,l);
            }
        }

        if(norm2_B_orth(k) >= (delta - mu(k,k-1)*mu(k,k-1)) * norm2_B_orth(k-1)){
            k++;
        }else{
            B.col(k-1).swap(B.col(k));
            U.col(k-1).swap(U.col(k));
            U_inv.row(k-1).swap(U_inv.row(k));

            Real mu_prime = mu(k,k-1);
            Real new_norm2_B = norm2_B_orth(k) + mu_prime*mu_prime * norm2_B_orth(k-1);
            mu(k,k-1) = mu_prime * norm2_B_orth(k-1) / new_norm2_B;
            norm2_B_orth(k) = norm2_B_orth(k) * norm2_B_orth(k-1) / new_norm2_B;
            norm2_B_orth(k-1) = new_norm2_B;

            for(int j = 0; j <= k-2; j++){
                std::swap(mu(k-1,j), mu(k,j));
            }

            for(int j = k+1; j < N; j++){
                Real t = mu(j,k);
                mu(j,k) = mu(j,k-1) - mu_prime * t;
                mu(j,k-1) = t + mu(k,k-1) * mu(j,k);
            }

            k = std::max(k-1, 1);
        }
    }

    prof.stop("LLL");
    return {U, U_inv};
}
#define X(Real) template  \
std::pair<Eigen::MatrixX<Integer>, Eigen::MatrixX<Integer>> LLL( \
    const Eigen::MatrixX<Real>&, \
    Real);                       \
REAL_SCALAR_TYPES
#undef X



template <typename Real>
std::vector<Eigen::VectorX<Integer>> EnumIntegerPoints(
        Eigen::MatrixX<Real> Q,
        Eigen::VectorX<Real> p,
        Real c
)
{
    auto& prof = util::Profiler::instance();
    prof.start("EnumIntegerPoints");

    int N = Q.rows();
    std::vector<Eigen::MatrixX<Real>> Q_list(N+1), Q_inv_list(N+1);
    std::vector<Eigen::MatrixX<Integer>> U_list(N+1), U_inv_list(N+1);
    std::vector<Eigen::MatrixX<Real>> U_real_list(N+1), U_inv_real_list(N+1);
    std::vector<Eigen::VectorX<Real>> A_invb_list(N+1);

    for(int i = N; i >= 1; i--) {
        if(i == N) Q_list[i] = Q;
        else Q_list[i] = Q_list[i+1].block(1,1, i,i);
        Q_inv_list[i] = math::linalg::InverseSPD(Q_list[i]);

        Eigen::MatrixX<Real> B = math::linalg::Cholesky(Q_inv_list[i]).transpose();
        auto [U_T, U_T_inv] = LLL(B, Real(0.75));
        Eigen::MatrixX<Integer> U = U_T.transpose();
        Eigen::MatrixX<Integer> U_inv = U_T_inv.transpose();
        Eigen::MatrixX<Real> U_real = U.unaryExpr([](Integer z) {
            return type::Int_to_Real<Real>(z);
        });
        Eigen::MatrixX<Real> U_inv_real = U_inv.unaryExpr([](Integer z) {
            return type::Int_to_Real<Real>(z);
        });

        Q_list[i]          = U_inv_real.transpose() * Q_list[i] * U_inv_real;
        Q_inv_list[i]      = U_real * Q_inv_list[i] * U_real.transpose();
        U_list[i]          = U;
        U_inv_list[i]      = U_inv;
        U_real_list[i]     = U_real;
        U_inv_real_list[i] = U_inv_real;
        A_invb_list[i]     = math::linalg::InverseSPD(Q_list[i].block(1,1, i-1,i-1).eval()) * Q_list[i].col(0).tail(i-1);
    }
    
    std::function<std::vector<Eigen::VectorX<Integer>>(int, Eigen::VectorX<Real>, Real)> _EnumIntegerPoints;
    _EnumIntegerPoints = [&](int N, Eigen::VectorX<Real> p, Real c)
    {
        std::vector<Eigen::VectorX<Integer>> ret;
        Eigen::MatrixX<Real> Q          = Q_list[N];
        Eigen::MatrixX<Real> Q_inv      = Q_inv_list[N];
        Eigen::MatrixX<Integer> U       = U_list[N];
        Eigen::MatrixX<Integer> U_inv   = U_inv_list[N];
        Eigen::MatrixX<Real> U_real     = U_real_list[N];
        Eigen::MatrixX<Real> U_inv_real = U_inv_real_list[N];
        Eigen::VectorX<Real> A_invb     = A_invb_list[N];

        if(N == 1){
            Real delta = math::sqrt(Q_inv(0,0) * c);
            Integer x0_min = math::ceil(p(0) - delta);
            Integer x0_max = math::floor(p(0) + delta);
            for(Integer x0 = x0_min; x0 <= x0_max; x0++){
                Eigen::Vector<Integer, 1> xvec;
                xvec(0) = x0;
                ret.push_back(xvec);
            }
        }else{
            p = U_real * p;
            // Q = U_inv_Real.transpose() * Q * U_inv_Real;
            // Q_inv = U_Real * Q_inv * U_Real.transpose();

            Real delta = math::sqrt(Q_inv(0,0) * c);
            Integer x0_min = math::ceil(p(0) - delta);
            Integer x0_max = math::floor(p(0) + delta);

            for(Integer x0 = x0_min; x0 <= x0_max; x0++){
                Eigen::VectorX<Integer> xvec(N);
                xvec(0) = x0;
                
                Real x0_real = type::Int_to_Real<Real>(x0);
                Eigen::VectorX<Real> b = (x0_real - p(0)) * Q.col(0).tail(N-1);
                Real next_c = c;
                Eigen::MatrixX<Real> next_p = p.tail(N-1);
                next_p -= (x0_real - p(0)) * A_invb;
                next_c -= Q(0,0) * (x0_real - p(0)) * (x0_real - p(0));
                next_c += b.dot((x0_real - p(0)) * A_invb);
                auto X_tail = _EnumIntegerPoints(N-1, next_p, next_c);
                for(auto xvec_tail : X_tail){
                    xvec.tail(N-1) = xvec_tail;

                    ret.push_back(U_inv * xvec);
                }
            }
        }
        return ret;
    };

    auto ret = _EnumIntegerPoints(N, p, c);
    
    prof.stop("EnumIntegerPoints");
    return ret;
}
#define X(Real) template \
std::vector<Eigen::VectorX<Integer>> EnumIntegerPoints( \
    Eigen::MatrixX<Real> Q, \
    Eigen::VectorX<Real> p, \
    Real c);
REAL_SCALAR_TYPES
#undef X




}