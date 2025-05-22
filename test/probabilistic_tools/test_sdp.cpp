#include <gtest/gtest.h>

#include <numeric>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>

#include "type.hpp"
#include "probabilistic_tools/sdp.hpp"

using namespace::su2compiler;

TEST(DiagBlockMatrix, Add) {
    int NumBlocks = 3;
    std::vector<int> BlockSizes = {3, 1, 2};
    std::vector<MatrixXR> mat1(NumBlocks), mat2(NumBlocks);
    for(int i = 0; i < NumBlocks; i++){
        mat1[i] = MatrixXR::Random(BlockSizes[i], BlockSizes[i]);
        mat2[i] = MatrixXR::Random(BlockSizes[i], BlockSizes[i]);
    }
    
    sdp::DiagBlockMatrix M1(mat1), M2(mat2);
    MatrixXR M1_Dense = M1.to_DenseMatrix();
    MatrixXR M2_Dense = M2.to_DenseMatrix();
    
    sdp::DiagBlockMatrix M1_add_M2 = M1 + M2;
    MatrixXR M1_add_M2_Dense = M1_add_M2.to_DenseMatrix();

    Real diff = (M1_add_M2_Dense - (M1_Dense + M2_Dense)).sum();
    EXPECT_NEAR(double(diff), 0.0, 1e-10);
}

TEST(DiagBlockMatrix, Sub) {
    int NumBlocks = 3;
    std::vector<int> BlockSizes = {3, 1, 2};
    std::vector<MatrixXR> mat1(NumBlocks), mat2(NumBlocks);
    for(int i = 0; i < NumBlocks; i++){
        mat1[i] = MatrixXR::Random(BlockSizes[i], BlockSizes[i]);
        mat2[i] = MatrixXR::Random(BlockSizes[i], BlockSizes[i]);
    }
    
    sdp::DiagBlockMatrix M1(mat1), M2(mat2);
    MatrixXR M1_Dense = M1.to_DenseMatrix();
    MatrixXR M2_Dense = M2.to_DenseMatrix();
    
    sdp::DiagBlockMatrix M1_sub_M2 = M1 - M2;
    MatrixXR M1_sub_M2_Dense = M1_sub_M2.to_DenseMatrix();

    Real diff = (M1_sub_M2_Dense - (M1_Dense - M2_Dense)).sum();
    EXPECT_NEAR(double(diff), 0.0, 1e-10);
}

TEST(DiagBlockMatrix, Mul) {
    int NumBlocks = 3;
    std::vector<int> BlockSizes = {3, 1, 2};
    std::vector<MatrixXR> mat1(NumBlocks), mat2(NumBlocks);
    for(int i = 0; i < NumBlocks; i++){
        mat1[i] = MatrixXR::Random(BlockSizes[i], BlockSizes[i]);
        mat2[i] = MatrixXR::Random(BlockSizes[i], BlockSizes[i]);
    }
    
    sdp::DiagBlockMatrix M1(mat1), M2(mat2);
    MatrixXR M1_Dense = M1.to_DenseMatrix();
    MatrixXR M2_Dense = M2.to_DenseMatrix();
    
    sdp::DiagBlockMatrix M1_mul_M2 = M1 * M2;
    MatrixXR M1_mul_M2_Dense = M1_mul_M2.to_DenseMatrix();

    Real diff = (M1_mul_M2_Dense - (M1_Dense * M2_Dense)).sum();
    EXPECT_NEAR(double(diff), 0.0, 1e-10);
}

TEST(DiagBlockMatrix, ScalarMul) {
    int NumBlocks = 3;
    std::vector<int> BlockSizes = {3, 1, 2};
    std::vector<MatrixXR> mat1(NumBlocks);
    for(int i = 0; i < NumBlocks; i++){
        mat1[i] = MatrixXR::Random(BlockSizes[i], BlockSizes[i]);
    }
    Real scalar = 100.0;
    
    sdp::DiagBlockMatrix M1(mat1);
    MatrixXR M1_Dense = M1.to_DenseMatrix();
    
    sdp::DiagBlockMatrix M1_mul_scalar = M1 * scalar;
    MatrixXR M1_mul_scalar_Dense = M1_mul_scalar.to_DenseMatrix();

    Real diff = (M1_mul_scalar_Dense - (M1_Dense * scalar)).sum();
    EXPECT_NEAR(double(diff), 0.0, 1e-10);
}

TEST(DiagBlockMatrix, inverse) {
    int NumBlocks = 3;
    std::vector<int> BlockSizes = {3, 1, 2};
    int NTotal = std::reduce(BlockSizes.begin(), BlockSizes.end());
    std::vector<MatrixXR> mat(NumBlocks);
    for(int i = 0; i < NumBlocks; i++){
        mat[i] = MatrixXR::Random(BlockSizes[i], BlockSizes[i]);
    }
    
    sdp::DiagBlockMatrix M(mat);
    sdp::DiagBlockMatrix M_inv = M.inverse();

    sdp::DiagBlockMatrix MM_inv = M * M_inv;
    Real diff = (MM_inv.to_DenseMatrix() - MatrixXR::Identity(NTotal, NTotal)).sum();
    EXPECT_NEAR(double(diff), 0.0, 1e-10);
}

TEST(DiagBlockMatrix, HSinner) {
    int NumBlocks = 3;
    std::vector<int> BlockSizes = {3, 1, 2};
    std::vector<MatrixXR> mat1(NumBlocks), mat2(NumBlocks);
    for(int i = 0; i < NumBlocks; i++){
        mat1[i] = MatrixXR::Random(BlockSizes[i], BlockSizes[i]);
        mat2[i] = MatrixXR::Random(BlockSizes[i], BlockSizes[i]);
    }
    
    sdp::DiagBlockMatrix M1(mat1), M2(mat2);
    MatrixXR M1_Dense = M1.to_DenseMatrix();
    MatrixXR M2_Dense = M2.to_DenseMatrix();
    
    Real value1 = sdp::HSinner(M1, M2);
    Real value2 = (M1_Dense * M2_Dense.transpose()).trace();

    Real diff = value1 - value2;
    EXPECT_NEAR(double(diff), 0.0, 1e-10);
}


TEST(SDP, example1) {
    int m = 3;
    int NumBlocks = 1;
    std::vector<int> BlockSizes = {2};

    VectorXR c(m);
    c << 48, -8, 20;

    std::vector<std::vector<MatrixXR>> Fmat(m+1);
    for(int i = 0; i <= m; i++){
        Fmat[i].resize(NumBlocks);
        for(int j = 0; j < NumBlocks; j++) Fmat[i][j] = MatrixXR::Zero(BlockSizes[j], BlockSizes[j]);
    }
    Fmat[0][0] << -11, 0,
                    0, 23;
    Fmat[1][0] << 10, 4,
                   4, 0;
    Fmat[2][0] << 0, 0,
                  0,-8;
    Fmat[3][0] << 0,-8,
                 -8,-2;
    std::vector<sdp::DiagBlockMatrix> F;
    for(int i = 0; i <= m; i++) F.push_back(sdp::DiagBlockMatrix(Fmat[i])); 

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


TEST(SDP, example2) {
    int m = 5;
    int NumBlocks = 4;
    std::vector<int> BlockSizes = {2,3,1,1};

    VectorXR c(m);
    c << 1.1, -10, 6.6 , 19 , 4.1;

    std::vector<std::vector<MatrixXR>> Fmat(m+1);
    for(int i = 0; i <= m; i++){
        Fmat[i].resize(NumBlocks);
        for(int j = 0; j < NumBlocks; j++) Fmat[i][j] = MatrixXR::Zero(BlockSizes[j], BlockSizes[j]);
    }
    
    Fmat[0][0] << -1.4, -3.2,
                  -3.2,-28;
    Fmat[0][1] << 15,  -12,    2.1,
                 -12,   16,   -3.8,
                  2.1, -3.8, 15;
    Fmat[0][2] << 1.8;
    Fmat[0][3] << -4.0;

    Fmat[1][0] << 0.5,  5.2,
                  5.2, -5.3;
    Fmat[1][1] << 7.8, -2.4,  6.0,
                 -2.4,  4.2,  6.5,
                 6.0,  6.5,  2.1;
    Fmat[1][2] << -4.5;
    Fmat[1][3] << -3.5;

    Fmat[2][0] << 1.7,  7.0,
                  7.0, -9.3 ;
    Fmat[2][1] << -1.9, -0.9, -1.3 ,
                  -0.9, -0.8, -2.1,
                  -1.3, -2.1,  4.0;
    Fmat[2][2] << -0.2;
    Fmat[2][3] << -3.7;

    Fmat[3][0] << 6.3, -7.5,
                 -7.5, -3.3;
    Fmat[3][1] << 0.2,  8.8,  5.4 ,
                  8.8,  3.4, -0.4,
                  5.4, -0.4,  7.5;
    Fmat[3][2] << -3.3;
    Fmat[3][3] << -4.0;

    Fmat[4][0] << -2.4, -2.5,
                  -2.5, -2.9;
    Fmat[4][1] << 3.4, -3.2, -4.5,
                 -3.2,  3.0, -4.8,
                 -4.5, -4.8,  3.6;
    Fmat[4][2] << 4.8;
    Fmat[4][3] << 9.7;

    Fmat[5][0] << -6.5, -5.4 ,
                  -5.4, -6.6;
    Fmat[5][1] << 6.7, -7.2, -3.6,
                 -7.2,  7.3, -3.0,
                 -3.6, -3.0, -1.4;
    Fmat[5][2] << 6.1;
    Fmat[5][3] << -1.5;


    std::vector<sdp::DiagBlockMatrix> F;
    for(int i = 0; i <= m; i++) F.push_back(sdp::DiagBlockMatrix(Fmat[i])); 

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