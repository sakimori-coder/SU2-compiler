#include <gtest/gtest.h>

#include <numeric>
#include <random>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>

#include "type.hpp"
#include "probabilistic_tools/sdp.hpp"

using namespace::su2compiler;

sdp::DenseMat DenseMatRandom(int m, int n) { return sdp::DenseMat::Random(m,n); }
sdp::DenseVec DenseVecRandom(int size) { return sdp::DenseVec::Random(size); }

sdp::SparseMat SparseMatRandom(int m, int n, double density = 0.5) {
    static std::mt19937_64 rng{std::random_device{}()};
    std::uniform_int_distribution<int>  rowDist(0, m-1);
    std::uniform_int_distribution<int>  colDist(0, n-1);
    std::uniform_real_distribution<double> valDist(-1.0, 1.0);

    const int nnz = static_cast<std::size_t>(m * n * density);

    std::vector<Eigen::Triplet<Real>> triplets;

    for(int i = 0; i < nnz; i++) {
        triplets.push_back({rowDist(rng), colDist(rng), valDist(rng)});
    }

    sdp::SparseMat A(m,n);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

sdp::SparseVec SparseVecRandom(int size, double density = 0.5) {
    static std::mt19937_64 rng{std::random_device{}()};
    std::uniform_int_distribution<int>  idxDist(0, size-1);
    std::uniform_real_distribution<double> valDist(-1.0, 1.0);

    const int nnz = static_cast<std::size_t>(size * density);

    sdp::SparseVec v(size);

    for(int i = 0; i < nnz; i++) {
        v.coeffRef(idxDist(rng)) = valDist(rng);
    }

    return v.pruned();
}


TEST(DiagBlockMatrix, transpose) {
    Real::set_default_prec(100);   // 30 digits precision

    int size = 6;
    std::vector<sdp::MatVecVar> Blocks;

    // DenseMat
    Blocks.push_back(DenseMatRandom(size, size));

    // SparseMat
    Blocks.push_back(SparseMatRandom(size, size));

    // DenseVec
    Blocks.push_back(DenseVecRandom(size));

    // SparseVec
    Blocks.push_back(SparseVecRandom(size));

    sdp::DiagBlockMatrix A(Blocks);
    sdp::DiagBlockMatrix A_transpose = A.transpose();
    MatrixXR A_transpose_Dense = A.to_DenseMatrix().transpose();

    EXPECT_TRUE(A_transpose.to_DenseMatrix().isApprox(A_transpose_Dense));

    std::vector<sdp::MatVecVar> Blocks_ret = A_transpose.get_Blocks();
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[0]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseMat>(Blocks_ret[1]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseVec>(Blocks_ret[2]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseVec>(Blocks_ret[3]));
}


TEST(DiagBlockMatrix, inverse) {
    Real::set_default_prec(100);   // 30 digits precision

    int size = 6;
    std::vector<sdp::MatVecVar> Blocks;

    sdp::SparseMat I(size, size);
    I.setIdentity();
    Real shift(1000);

    // DenseMat
    Blocks.push_back(sdp::DenseMat(DenseMatRandom(size, size) + shift * I));

    // SparseMat
    Blocks.push_back(sdp::SparseMat(SparseMatRandom(size, size) + shift * I));

    // DenseVec
    Blocks.push_back(sdp::DenseVec(DenseVecRandom(size) + shift * VectorXR::Ones(size)));

    sdp::DiagBlockMatrix A(Blocks);
    sdp::DiagBlockMatrix A_inv = A.inverse();
    MatrixXR A_inv_Dense = A.to_DenseMatrix().inverse();

    EXPECT_TRUE(A_inv.to_DenseMatrix().isApprox(A_inv_Dense));

    std::vector<sdp::MatVecVar> Blocks_ret = A_inv.get_Blocks();
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[0]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[1]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseVec>(Blocks_ret[2]));

    std::cout << A_inv.to_DenseMatrix() * A.to_DenseMatrix() << std::endl;
}


TEST(DiagBlockMatrix, Add) {
    Real::set_default_prec(100);   // 30 digits precision

    int size = 6;
    std::vector<sdp::MatVecVar> Blocks1, Blocks2;

    // DenseMat + DenseMat
    Blocks1.push_back(DenseMatRandom(size, size));
    Blocks2.push_back(DenseMatRandom(size, size));
    
    // SparseMat + DenseMat
    Blocks1.push_back(SparseMatRandom(size, size));
    Blocks2.push_back(DenseMatRandom(size, size));

    // DenseMat + SparseMat
    Blocks1.push_back(DenseMatRandom(size, size));
    Blocks2.push_back(SparseMatRandom(size, size));

    // SparseMat + SparseMat
    Blocks1.push_back(SparseMatRandom(size, size));
    Blocks2.push_back(SparseMatRandom(size, size));

    // DenseVec + DenseVec
    Blocks1.push_back(DenseVecRandom(size));
    Blocks2.push_back(DenseVecRandom(size));

    // SparseVec + DenseVec
    Blocks1.push_back(SparseVecRandom(size));
    Blocks2.push_back(DenseVecRandom(size));

    // DenseVec + SparseVec
    Blocks1.push_back(DenseVecRandom(size));
    Blocks2.push_back(SparseVecRandom(size));

    // SparseVec + SparseVec
    Blocks1.push_back(SparseVecRandom(size));
    Blocks2.push_back(SparseVecRandom(size));


    sdp::DiagBlockMatrix A(Blocks1), B(Blocks2);
    sdp::DiagBlockMatrix C = A + B;
    MatrixXR A_Dense = A.to_DenseMatrix();
    MatrixXR B_Dense = B.to_DenseMatrix();
    MatrixXR C_Dense = A_Dense + B_Dense;

    EXPECT_TRUE(C.to_DenseMatrix().isApprox(C_Dense, Real(1e-20)));

    std::vector<sdp::MatVecVar> Blocks_ret = C.get_Blocks();
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[0]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[1]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[2]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseMat>(Blocks_ret[3]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseVec>(Blocks_ret[4]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseVec>(Blocks_ret[5]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseVec>(Blocks_ret[6]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseVec>(Blocks_ret[7]));
}


TEST(DiagBlockMatrix, Sub) {
    Real::set_default_prec(100);   // 30 digits precision

    int size = 6;
    std::vector<sdp::MatVecVar> Blocks1, Blocks2;

    // DenseMat - DenseMat
    Blocks1.push_back(DenseMatRandom(size, size));
    Blocks2.push_back(DenseMatRandom(size, size));
    
    // SparseMat - DenseMat
    Blocks1.push_back(SparseMatRandom(size, size));
    Blocks2.push_back(DenseMatRandom(size, size));

    // DenseMat - SparseMat
    Blocks1.push_back(DenseMatRandom(size, size));
    Blocks2.push_back(SparseMatRandom(size, size));

    // SparseMat - SparseMat
    Blocks1.push_back(SparseMatRandom(size, size));
    Blocks2.push_back(SparseMatRandom(size, size));

    // DenseVec - DenseVec
    Blocks1.push_back(DenseVecRandom(size));
    Blocks2.push_back(DenseVecRandom(size));

    // SparseVec - DenseVec
    Blocks1.push_back(SparseVecRandom(size));
    Blocks2.push_back(DenseVecRandom(size));

    // DenseVec - SparseVec
    Blocks1.push_back(DenseVecRandom(size));
    Blocks2.push_back(SparseVecRandom(size));

    // SparseVec - SparseVec
    Blocks1.push_back(SparseVecRandom(size));
    Blocks2.push_back(SparseVecRandom(size));

    sdp::DiagBlockMatrix A(Blocks1), B(Blocks2);
    sdp::DiagBlockMatrix C = A - B;
    MatrixXR A_Dense = A.to_DenseMatrix();
    MatrixXR B_Dense = B.to_DenseMatrix();
    MatrixXR C_Dense = A_Dense - B_Dense;

    EXPECT_TRUE(C.to_DenseMatrix().isApprox(C_Dense, Real(1e-20)));

    std::vector<sdp::MatVecVar> Blocks_ret = C.get_Blocks();
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[0]));
    // EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[1]));
    // EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[2]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseMat>(Blocks_ret[1]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseMat>(Blocks_ret[2]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseMat>(Blocks_ret[3]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseVec>(Blocks_ret[4]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseVec>(Blocks_ret[5]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseVec>(Blocks_ret[6]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseVec>(Blocks_ret[7]));
}


TEST(DiagBlockMatrix, Mul) {
    Real::set_default_prec(100);   // 30 digits precision

    int size = 6;
    std::vector<sdp::MatVecVar> Blocks1, Blocks2;

    // DenseMat * DenseMat
    Blocks1.push_back(DenseMatRandom(size, size));
    Blocks2.push_back(DenseMatRandom(size, size));
    
    // SparseMat * DenseMat
    Blocks1.push_back(SparseMatRandom(size, size));
    Blocks2.push_back(DenseMatRandom(size, size));

    // DenseMat * SparseMat
    Blocks1.push_back(DenseMatRandom(size, size));
    Blocks2.push_back(SparseMatRandom(size, size));

    // SparseMat * SparseMat
    Blocks1.push_back(SparseMatRandom(size, size));
    Blocks2.push_back(SparseMatRandom(size, size));

    // DenseVec * DenseVec
    Blocks1.push_back(DenseVecRandom(size));
    Blocks2.push_back(DenseVecRandom(size));

    // SparseVec * DenseVec
    Blocks1.push_back(SparseVecRandom(size));
    Blocks2.push_back(DenseVecRandom(size));

    // DenseVec * SparseVec
    Blocks1.push_back(DenseVecRandom(size));
    Blocks2.push_back(SparseVecRandom(size));

    // SparseVec * SparseVec
    Blocks1.push_back(SparseVecRandom(size));
    Blocks2.push_back(SparseVecRandom(size));

    sdp::DiagBlockMatrix A(Blocks1), B(Blocks2);
    sdp::DiagBlockMatrix C = A * B;
    MatrixXR A_Dense = A.to_DenseMatrix();
    MatrixXR B_Dense = B.to_DenseMatrix();
    MatrixXR C_Dense = A_Dense * B_Dense;

    EXPECT_TRUE(C.to_DenseMatrix().isApprox(C_Dense, Real(1e-20)));

    std::vector<sdp::MatVecVar> Blocks_ret = C.get_Blocks();
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[0]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[1]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[2]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseMat>(Blocks_ret[3]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseVec>(Blocks_ret[4]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseVec>(Blocks_ret[5]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseVec>(Blocks_ret[6]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseVec>(Blocks_ret[7]));
}


TEST(DiagBlockMatrix, ScalarMul) {
    Real::set_default_prec(100);   // 30 digits precision

    int size = 6;
    std::vector<sdp::MatVecVar> Blocks;

    // scalar * DenseMat
    Blocks.push_back(DenseMatRandom(size, size));

    // scalar * SparseMat
    Blocks.push_back(SparseMatRandom(size, size));

    // scalar * DenseVec
    Blocks.push_back(DenseVecRandom(size));

    // scalar * SparseVec
    Blocks.push_back(SparseVecRandom(size));

    Real c = 100;
    sdp::DiagBlockMatrix A(Blocks);
    sdp::DiagBlockMatrix cA = c * A;
    MatrixXR A_Dense = A.to_DenseMatrix();
    MatrixXR cA_Dense = c * A_Dense;

    EXPECT_TRUE(cA.to_DenseMatrix().isApprox(cA_Dense, Real(1e-20)));

    std::vector<sdp::MatVecVar> Blocks_ret = cA.get_Blocks();
    EXPECT_TRUE(std::holds_alternative<sdp::DenseMat>(Blocks_ret[0]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseMat>(Blocks_ret[1]));
    EXPECT_TRUE(std::holds_alternative<sdp::DenseVec>(Blocks_ret[2]));
    EXPECT_TRUE(std::holds_alternative<sdp::SparseVec>(Blocks_ret[3]));
}




TEST(DiagBlockMatrix, HSinner) {
    Real::set_default_prec(100);   // 30 digits precision

    int size = 6;
    std::vector<sdp::MatVecVar> Blocks1, Blocks2;

    // DenseMat * DenseMat
    Blocks1.push_back(DenseMatRandom(size, size));
    Blocks2.push_back(DenseMatRandom(size, size));
    
    // SparseMat * DenseMat
    Blocks1.push_back(SparseMatRandom(size, size));
    Blocks2.push_back(DenseMatRandom(size, size));

    // DenseMat * SparseMat
    Blocks1.push_back(DenseMatRandom(size, size));
    Blocks2.push_back(SparseMatRandom(size, size));

    // SparseMat * SparseMat
    Blocks1.push_back(SparseMatRandom(size, size));
    Blocks2.push_back(SparseMatRandom(size, size));

    // DenseVec * DenseVec
    Blocks1.push_back(DenseVecRandom(size));
    Blocks2.push_back(DenseVecRandom(size));

    // SparseVec * DenseVec
    Blocks1.push_back(SparseVecRandom(size));
    Blocks2.push_back(DenseVecRandom(size));

    // DenseVec * SparseVec
    Blocks1.push_back(DenseVecRandom(size));
    Blocks2.push_back(SparseVecRandom(size));

    // SparseVec * SparseVec
    Blocks1.push_back(SparseVecRandom(size));
    Blocks2.push_back(SparseVecRandom(size));

    sdp::DiagBlockMatrix A(Blocks1), B(Blocks2);
    Real val1 = sdp::HSinner(A, B);
    MatrixXR A_Dense = A.to_DenseMatrix();
    MatrixXR B_Dense = B.to_DenseMatrix();
    Real val2 = A_Dense.cwiseProduct(B_Dense).sum();

    EXPECT_NEAR(double(val1), double(val2), 1e-13);
}


TEST(SDP, example1) {
    int m = 3;
    int NumBlocks = 1;
    std::vector<int> BlockSizes = {2};

    VectorXR c(m);
    c << 48, -8, 20;

    std::vector<std::vector<MatrixXR>> F_Blocks(m+1);
    for(int i = 0; i <= m; i++){
        F_Blocks[i].resize(NumBlocks);
        for(int j = 0; j < NumBlocks; j++) F_Blocks[i][j] = MatrixXR::Zero(BlockSizes[j], BlockSizes[j]);
    }
    F_Blocks[0][0] << -11, 0,
                       0, 23;
    F_Blocks[1][0] << 10, 4,
                       4, 0;
    F_Blocks[2][0] << 0, 0,
                      0,-8;
    F_Blocks[3][0] << 0,-8,
                     -8,-2;

    std::vector<std::vector<sdp::MatVecVar>> F_BlocksVar(m+1);
    for(int i = 0; i < m+1; i++) {
        for(auto& blk : F_Blocks[i]) F_BlocksVar[i].push_back(blk);
    }
    std::vector<sdp::DiagBlockMatrix> F;
    for(int i = 0; i <= m; i++) F.push_back(sdp::DiagBlockMatrix(F_BlocksVar[i])); 

    // for(int i = 0; i <= m; i++) std::cout << F[i] << std::endl;

    Real lambda = 1e4;
    Real betaStar = 0.1;
    Real betaBar = 0.3;
    Real gamma = 0.9;
    Real epsilon1 = 1e-30;
    Real epsilon2 = 1e-30;
    int maxITERATION = 200;
    bool OUTPUT_HISTORY = true;

    VectorXR x = sdp::SDP(
        c,
        F,
        lambda,
        betaStar,
        betaBar,
        gamma,
        epsilon1,
        epsilon2,
        maxITERATION,
        OUTPUT_HISTORY
    );

    std::cout << c.dot(x) << std::endl;
}


// TEST(SDP, example2) {
//     mpfr::mpreal::set_default_prec(1024);
//     int m = 5;
//     int NumBlocks = 3;
//     std::vector<int> BlockSizes = {2,3,-2};

//     VectorXR c(m);
//     c << 1.1, -10, 6.6 , 19 , 4.1;

//     std::vector<std::vector<MatrixXR>> F_Blocks(m+1);
//     for(int i = 0; i <= m; i++){
//         F_Blocks[i].resize(NumBlocks);
//         for(int j = 0; j < NumBlocks; j++) {
//             if(BlockSizes[j] > 0) F_Blocks[i][j] = MatrixXR::Zero(BlockSizes[j], BlockSizes[j]);
//             else                  F_Blocks[i][j] = VectorXR::Zero(-BlockSizes[j]);
//         }
//     }
    
//     F_Blocks[0][0] << -1.4, -3.2,
//                   -3.2,-28;
//     F_Blocks[0][1] << 15,  -12,    2.1,
//                  -12,   16,   -3.8,
//                   2.1, -3.8, 15;
//     F_Blocks[0][2] << 1.8, -4.0;

//     F_Blocks[1][0] << 0.5,  5.2,
//                   5.2, -5.3;
//     F_Blocks[1][1] << 7.8, -2.4,  6.0,
//                  -2.4,  4.2,  6.5,
//                  6.0,  6.5,  2.1;
//     F_Blocks[1][2] << -4.5, -3.5;

//     F_Blocks[2][0] << 1.7,  7.0,
//                   7.0, -9.3 ;
//     F_Blocks[2][1] << -1.9, -0.9, -1.3 ,
//                   -0.9, -0.8, -2.1,
//                   -1.3, -2.1,  4.0;
//     F_Blocks[2][2] << -0.2, -3.7;

//     F_Blocks[3][0] << 6.3, -7.5,
//                  -7.5, -3.3;
//     F_Blocks[3][1] << 0.2,  8.8,  5.4 ,
//                   8.8,  3.4, -0.4,
//                   5.4, -0.4,  7.5;
//     F_Blocks[3][2] << -3.3, -4.0;

//     F_Blocks[4][0] << -2.4, -2.5,
//                   -2.5, -2.9;
//     F_Blocks[4][1] << 3.4, -3.2, -4.5,
//                  -3.2,  3.0, -4.8,
//                  -4.5, -4.8,  3.6;
//     F_Blocks[4][2] << 4.8, 9.7;

//     F_Blocks[5][0] << -6.5, -5.4 ,
//                   -5.4, -6.6;
//     F_Blocks[5][1] << 6.7, -7.2, -3.6,
//                  -7.2,  7.3, -3.0,
//                  -3.6, -3.0, -1.4;
//     F_Blocks[5][2] << 6.1, -1.5;


//     std::vector<sdp::DiagBlockMatrix> F;
//     for(int i = 0; i <= m; i++) F.push_back(sdp::DiagBlockMatrix(F_Blocks[i])); 

//     // for(int i = 0; i <= m; i++) std::cout << F[i] << std::endl;

//     Real lambda = 1e4;
//     Real betaStar = 0.1;
//     Real betaBar = 0.3;
//     Real gamma = 0.9;
//     Real epsilon1 = 1e-30;
//     Real epsilon2 = 1e-30;
//     int maxITERATION = 200;
//     bool OUTPUT_HISTORY = true;

//     VectorXR x = sdp::SDP(
//         c,
//         F,
//         lambda,
//         betaStar,
//         betaBar,
//         gamma,
//         epsilon1,
//         epsilon2,
//         maxITERATION,
//         OUTPUT_HISTORY
//     );

//     std::cout << c.dot(x) << std::endl;
// }