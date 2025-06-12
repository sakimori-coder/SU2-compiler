#include <gtest/gtest.h>

#include "core/type.hpp"
#include "math/sdp.hpp"

using namespace su2compiler;
using namespace su2compiler::math::sdp;

TEST(SdpTest, example1)
{
    Real::set_default_prec(256);

    int M = 3;
    int NumBlocks = 1;
    std::vector<int> BlockSizes = {2};

    VectorR c(M);
    c << 48, -8, 20;

    std::vector<std::vector<MatrixR>> F_Blocks(M+1);
    for(int i = 0; i <= M; i++){
        F_Blocks[i].resize(NumBlocks);
        for(int j = 0; j < NumBlocks; j++) F_Blocks[i][j] = MatrixR::Zero(BlockSizes[j], BlockSizes[j]);
    }
    F_Blocks[0][0] << -11, 0,
                       0, 23;
    F_Blocks[1][0] << 10, 4,
                       4, 0;
    F_Blocks[2][0] << 0, 0,
                      0,-8;
    F_Blocks[3][0] << 0,-8,
                     -8,-2;

    std::vector<SparseDiagBlock> F(M+1);
    for(int i = 0; i <= M; i++) {
        std::vector<SparseMatrixR> Fi_Blocks(NumBlocks);
        for(int j = 0; j < NumBlocks; j++) {
            Fi_Blocks[j] = F_Blocks[i][j].sparseView();
        } 
        F[i] = SparseDiagBlock(Fi_Blocks);
    }

    Options opt;
    opt.MaxIteration = 200;
    opt.lambda       = 1e4;
    opt.betaStar     = 0.1;
    opt.betaBar      = 0.3;
    opt.gamma        = 0.9;
    opt.epsilon1     = 1e-30;
    opt.epsilon2     = 1e-30;
    opt.OUTPUT_HISTORY = true;
    

    Results ret = solve(
        c,
        F,
        opt
    );

    EXPECT_TRUE(ret.status == STATUS::SUCCESS);
}



TEST(SdpTest, example2) {
    Real::set_default_prec(256);
    
    int M = 5;
    int NumBlocks = 3;
    std::vector<int> BlockSizes = {2,3,-2};

    VectorR c(M);
    c << 1.1, -10, 6.6 , 19 , 4.1;

    std::vector<std::vector<MatrixR>> F_Blocks(M+1);
    for(int i = 0; i <= M; i++){
        F_Blocks[i].resize(NumBlocks);
        for(int j = 0; j < NumBlocks; j++) {
            if(BlockSizes[j] > 0) F_Blocks[i][j] = MatrixR::Zero(BlockSizes[j], BlockSizes[j]);
            else                  F_Blocks[i][j] = VectorR::Zero(-BlockSizes[j]);
        }
    }
    
    F_Blocks[0][0] << -1.4, -3.2,
                    -3.2,-28;
    F_Blocks[0][1] << 15,  -12,    2.1,
                 -12,   16,   -3.8,
                  2.1, -3.8, 15;
    F_Blocks[0][2] << 1.8, -4.0;

    F_Blocks[1][0] << 0.5,  5.2,
                  5.2, -5.3;
    F_Blocks[1][1] << 7.8, -2.4,  6.0,
                 -2.4,  4.2,  6.5,
                 6.0,  6.5,  2.1;
    F_Blocks[1][2] << -4.5, -3.5;

    F_Blocks[2][0] << 1.7,  7.0,
                  7.0, -9.3 ;
    F_Blocks[2][1] << -1.9, -0.9, -1.3 ,
                  -0.9, -0.8, -2.1,
                  -1.3, -2.1,  4.0;
    F_Blocks[2][2] << -0.2, -3.7;

    F_Blocks[3][0] << 6.3, -7.5,
                 -7.5, -3.3;
    F_Blocks[3][1] << 0.2,  8.8,  5.4 ,
                  8.8,  3.4, -0.4,
                  5.4, -0.4,  7.5;
    F_Blocks[3][2] << -3.3, -4.0;

    F_Blocks[4][0] << -2.4, -2.5,
                  -2.5, -2.9;
    F_Blocks[4][1] << 3.4, -3.2, -4.5,
                 -3.2,  3.0, -4.8,
                 -4.5, -4.8,  3.6;
    F_Blocks[4][2] << 4.8, 9.7;

    F_Blocks[5][0] << -6.5, -5.4 ,
                  -5.4, -6.6;
    F_Blocks[5][1] << 6.7, -7.2, -3.6,
                 -7.2,  7.3, -3.0,
                 -3.6, -3.0, -1.4;
    F_Blocks[5][2] << 6.1, -1.5;


    std::vector<SparseDiagBlock> F(M+1);
    for(int i = 0; i <= M; i++) {
        std::vector<SparseMatrixR> Fi_Blocks(NumBlocks);
        for(int j = 0; j < NumBlocks; j++) {
            Fi_Blocks[j] = F_Blocks[i][j].sparseView();
        } 
        F[i] = SparseDiagBlock(Fi_Blocks);
    }

    Options opt;
    opt.MaxIteration = 200;
    opt.lambda       = 1e4;
    opt.betaStar     = 0.1;
    opt.betaBar      = 0.3;
    opt.gamma        = 0.9;
    opt.epsilon1     = 1e-30;
    opt.epsilon2     = 1e-30;
    opt.OUTPUT_HISTORY = true;
    

    Results ret = solve(
        c,
        F,
        opt
    );

    EXPECT_TRUE(ret.status == STATUS::SUCCESS);

}