#include "core/mix_su2.hpp"

#include <vector>
#include <Eigen/Core>

#include "core/type.hpp"
#include "core/su2.hpp"


namespace su2compiler
{

using Eigen::Matrix4;
using Eigen::VectorX;

// template <typename Real>
// void MixSU2<Real>::compute_optimal_prob(
//     const SU2<Real>& targetV,
//     Real lambda,
//     Real betaStar,
//     Real betaBar,
//     Real gamma,
//     Real epsilon1,
//     Real epsilon2,
//     int  MaxIteration,
//     bool PUTPUT_HISTORY = false)
// {
//     const int N = U.size();
//     const int M = 10 + N;        // number of variables in SDP
//     const int L = 4 + 4 + N + 1; // dimension of F

//     std::vector<sdp_tool::BlockMatrix<Real>> F;

//     sdp_tool::BlockMatrix<Real> I(
//         Matrix4<Real>::Identity(),
//         Matrix4<Real>::Identity(),
//         VectorX<Real>::Ones(),
//         Real(1)
//     );

//     sdp_tool::BlockMatrix<Real> X = I * lambda;
//     sdp_tool::BlockMatrix<Real> Y = I * lambda;
//     VectorX<Real> x = VectorX<Real>(M);
//     Real mu = HSinner(X, Y) / L;
//     int Iteration = 0;

//     while(Iteration <= MaxIteration) {
//         sdp_tool::BlockMatrix<Real> Xinv;

//         Eigen::MatrixX<Real> B(M,M);
//         for(int i = 0; i < M; i++) {
//             sdp_tool::BlockMatrix<Real> left = Xinv * F[i+1] * Y;
//             for( int j = 0; j <= i; j++) {
//                 if(i < 10) {

//                 } else if(j < 10) {
                    
//                 } else {
//                     B(i,j) = B(j,i) = left.blk2.cwiseProduct(F[j].blk2) - left.blk4;
//                 }
//             }
//         }
//     }
    
    
    
// }


}