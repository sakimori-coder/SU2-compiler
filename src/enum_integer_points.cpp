#include "enum_integer_points.hpp"

#include <tuple>
#include <utility>
#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <chrono>

#include "type.hpp"

namespace su2compiler{


MatrixXR GSO(const MatrixXR& B)
{
    int M = B.rows();
    int N = B.cols();
    MatrixXR B_orth(M, N);
    auto mu = [&](int i, int j) {
        return B.col(i).dot(B_orth.col(j)) / B_orth.col(j).squaredNorm();
    }; 

    for(int i = 0; i < N; i++){
        B_orth.col(i) = B.col(i);
        for(int j = 0; j < i; j++){
            B_orth.col(i) -= mu(i, j) * B_orth.col(j);
        }
    }
    return B_orth;
}



std::tuple<MatrixXI, MatrixXI> LLL(MatrixXR B, Real delta)
{
    int M = B.rows();
    int N = B.cols();
    MatrixXR B_orth(M, N);
    MatrixXR mu(N, N);
    VectorXR norm2_B_orth(N);

    B_orth = GSO(B);
    for(int i = 0; i < N; i++) norm2_B_orth[i] = B_orth.col(i).squaredNorm();
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            mu(i, j) = B.col(i).dot(B_orth.col(j)) / B_orth.col(j).squaredNorm();
        }
    }

    MatrixXI U = MatrixXI::Identity(N, N);
    MatrixXI U_inv = MatrixXI::Identity(N, N);

    int k = 1;
    while(k < N){
        // std::cout << k << std::endl;
        for(int j = k-1; j >= 0; j--){
            if(abs(mu(k,j)) > 0.5){
                Integer q = round_mpreal(mu(k,j));
                Real q_real = Real(q.get_mpz_t());
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
    return {U, U_inv};
}



MatrixXR cholesky(const MatrixXR& A){
    int N = A.rows();
    MatrixXR L = MatrixXR::Zero(N,N);

    for(int i = 0; i < N; i++){
        for(int j = 0; j <= i; j++){
            Real sum(0);
            for(int k = 0; k < j; k++){
                sum += L(i,k) * L(j,k);
            }

            if(i == j){
                Real val = A(i,i) - sum;
                L(i,j) = sqrt(val);
            }else{
                L(i,j) = (A(i,j) - sum) / L(j,j);
            }
        }
    }
    
    return L;
}



std::vector<VectorXI> EnumIntegerPoints(MatrixXR Q, VectorXR p, Real c)
{
    int N = Q.rows();
    std::vector<MatrixXR> Q_list(N+1), Q_inv_list(N+1);
    std::vector<MatrixXI> U_list(N+1), U_inv_list(N+1);
    std::vector<VectorXR> A_invb_list(N+1);

    std::chrono::system_clock::time_point start, end;
    double time;
    start = std::chrono::system_clock::now();
    for(int i = N; i >= 1; i--){
        if(i == N) Q_list[i] = Q;
        else Q_list[i] = Q_list[i+1].block(1,1, i,i);
        Q_inv_list[i] = Q_list[i].inverse();

        MatrixXR B = cholesky(Q_inv_list[i]).transpose();
        auto [U_T, U_T_inv] = LLL(B, Real(0.75));
        // if(i != 8) U_T = U_T_inv = MatrixXI::Identity(i,i);
        MatrixXI U = U_T.transpose();
        MatrixXI U_inv = U_T_inv.transpose();
        MatrixXR U_Real = U.cast<Real>();
        MatrixXR U_inv_Real = U_inv.cast<Real>();
        Q_list[i] = U_inv_Real.transpose() * Q_list[i] * U_inv_Real;
        Q_inv_list[i] = U_Real * Q_inv_list[i] * U_Real.transpose();
        U_list[i] = U;
        U_inv_list[i] = U_inv;
        A_invb_list[i] = Q_list[i].block(1,1, i-1,i-1).inverse() * Q_list[i].col(0).tail(i-1);
    }
    end = std::chrono::system_clock::now();
    time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
    // std::cout << "LLL : " << time << "[ms]" << std::endl;

    std::function<std::vector<VectorXI>(int, VectorXR, Real)> _EnumIntegerPoints;
    _EnumIntegerPoints = [&_EnumIntegerPoints, &Q_list, &Q_inv_list, &U_list, &U_inv_list, &A_invb_list]
    (int N, VectorXR p, Real c)
    {
        std::vector<VectorXI> ret;
        MatrixXR Q = Q_list[N];
        MatrixXR Q_inv = Q_inv_list[N];
        MatrixXI U = U_list[N];
        MatrixXI U_inv = U_inv_list[N];
        MatrixXR U_Real = U.cast<Real>();
        MatrixXR U_inv_Real = U_inv.cast<Real>();
        VectorXR A_invb = A_invb_list[N];
        // std::cout << N << std::endl;

        if(N == 1){
            Real delta = sqrt(Q_inv(0,0) * c);
            Integer x0_min = ceil_mpreal(p(0) - delta);
            Integer x0_max = floor_mpreal(p(0) + delta);
            for(Integer x0 = x0_min; x0 <= x0_max; x0++){
                Eigen::Vector<Integer, 1> xvec;
                xvec(0) = x0;
                ret.push_back(xvec);
            }
        }else{
            p = U_Real * p;
            // Q = U_inv_Real.transpose() * Q * U_inv_Real;
            // Q_inv = U_Real * Q_inv * U_Real.transpose();

            Real delta = sqrt(Q_inv(0,0) * c);
            Integer x0_min = ceil_mpreal(p(0) - delta);
            Integer x0_max = floor_mpreal(p(0) + delta);

            for(Integer x0 = x0_min; x0 <= x0_max; x0++){
                VectorXI xvec(N);
                xvec(0) = x0;
                
                Real x0_real = Real(x0.get_mpz_t());
                VectorXR b = (x0_real - p(0)) * Q.col(0).tail(N-1);
                Real next_c = c;
                MatrixXR next_p = p.tail(N-1);
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

    start = std::chrono::system_clock::now();
    auto ret = _EnumIntegerPoints(N, p, c);
    end = std::chrono::system_clock::now();
    time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
    // std::cout << "Enum : " << time << "[ms]" << std::endl;
    
    return ret;
}


}