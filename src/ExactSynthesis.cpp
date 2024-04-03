#include "ExactSynthesis.hpp"
#include "U2_ZOmega.hpp"

namespace SU2_Compiler
{
    std::pair<std::string, std::vector<std::vector<int>>> 
    Clifford[24] = {{"", 
                     {{1, 0, 0},
                      {0, 1, 0},
                      {0, 0, 1}}},
                    
                    {"S",
                     {{0, -1, 0},
                      {1, 0, 0},
                      {0, 0, 1}}},
                    
                    {"H",
                     {{0, 0, 1},
                      {0, -1, 0},
                      {1, 0, 0}}},
                    
                    {"SS",
                     {{-1, 0, 0},
                      {0, -1, 0},
                      {0, 0, 1}}},
                    
                    {"HS",
                     {{0, 0, 1},
                      {-1, 0, 0},
                      {0, -1, 0}}},
                    
                    {"SH",
                     {{0, 1, 0},
                     {0, 0, 1},
                     {1, 0, 0}}},
                    
                    {"SSS",
                     {{0, 1, 0},
                      {-1, 0, 0},
                      {0, 0, 1}}},
                    
                    {"HSS",
                     {{0, 0, 1},
                      {0, 1, 0},
                      {-1, 0, 0}}},
                    
                    {"SHS",
                     {{1, 0, 0},
                     {0, 0, 1},
                     {0, -1, 0}}},
                    
                    {"SSH",
                     {{0, 0, -1},
                      {0, 1, 0},
                      {1, 0, 0}}},
                    
                    {"HSH",
                     {{1, 0, 0},
                      {0, 0, -1},
                      {0, 1, 0}}},
                    
                    {"HSSS",
                     {{0, 0, 1},
                      {1, 0, 0},
                      {0, 1, 0}}},
                    
                    {"SHSS",
                     {{0, -1, 0},
                      {0, 0, 1},
                      {-1, 0, 0}}},
                    
                    {"SSHS",
                     {{0, 0, -1},
                      {1, 0, 0},
                      {0, -1, 0}}},
                    
                    {"HSHS",
                     {{0, -1, 0},
                      {0, 0, -1},
                      {1, 0, 0}}},
                    
                    {"HSSH",
                     {{1, 0, 0},
                     {0, -1, 0},
                     {0, 0, -1}}},
                    
                    {"SHSSS",
                     {{-1, 0, 0},
                      {0, 0, 1},
                      {0, 1, 0}}},
                    
                    {"SSHSS",
                     {{0, 0, -1},
                      {0, -1, 0},
                      {-1, 0, 0}}},
                    
                    {"HSHSS",
                     {{-1, 0, 0},
                      {0, 0, -1},
                      {0, -1, 0}}},
                    
                    {"HSSHS",
                     {{0, -1, 0},
                      {-1, 0, 0},
                      {0, 0, -1}}},
                    
                    {"SHSSH",
                     {{0, 1, 0},
                      {1, 0, 0},
                      {0, 0, -1}}},
                    
                    {"SSHSSS",
                     {{0, 0, -1},
                      {-1, 0, 0},
                      {0, 1, 0}}},
                    
                    {"HSHSSS",
                     {{0, 1, 0},
                      {0, 0, -1},
                      {-1, 0, 0}}},
                    
                    {"HSSHSS",
                     {{-1, 0, 0},
                      {0, 1, 0},
                      {0, 0, -1}}}
                   };


    void reduction(Mat3ZRoot2& A, int& k)
    {
        while(true){
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++){
                    ITYPE a = A(i,j).a;
                    ITYPE b = A(i,j).b;
                    if(a % 2) return;
                }
            }

            k--;
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++){
                ITYPE a = A(i,j).a;
                ITYPE b = A(i,j).b;

                A(i,j).a = b;
                A(i,j).b = a / 2;
                }
            }
        }
    }


    std::pair<Mat3ZRoot2, int> U2_to_SO3(const Mat2ZOmega U, int k)
    {
        Mat2ZOmega U_dag;
        U_dag << conj(U(0,0)), conj(U(1,0)),
                 conj(U(0,1)), conj(U(1,1));
        
        Mat3ZRoot2 U_SO3;
        int k_SO3 = 2 * k + 1;

        ZOmega zero(0,0,0,0);
        ZOmega one(0,0,0,1);
        ZOmega ii(0,1,0,0);
        ZOmega sqrt2(-1,0,1,0);
        
        Eigen::Matrix<ZOmega, 2, 2> X, Y, Z;
        X << zero, one, 
            one , zero;
        Y << zero, -ii,
            ii, zero;
        Z << one , zero,
            zero, -one;
        

        Mat2ZOmega UXU_dag, UYU_dag, UZU_dag;
        UXU_dag = sqrt2 * U * X * U_dag;
        UYU_dag = sqrt2 * U * Y * U_dag;
        UZU_dag = sqrt2 * U * Z * U_dag;

        U_SO3(0,0) = UXU_dag(0,1).real();
        U_SO3(0,1) = UYU_dag(0,1).real();
        U_SO3(0,2) = UZU_dag(0,1).real();

        U_SO3(1,0) = -UXU_dag(0,1).imag();
        U_SO3(1,1) = -UYU_dag(0,1).imag();
        U_SO3(1,2) = -UZU_dag(0,1).imag();

        U_SO3(2,0) = UXU_dag(0,0).real();
        U_SO3(2,1) = UYU_dag(0,0).real();
        U_SO3(2,2) = UZU_dag(0,0).real();

        reduction(U_SO3, k_SO3);

        return {U_SO3, k_SO3};
    }
    
    Eigen::Matrix3i parity(const Mat3ZRoot2& A)
    {
        Eigen::Matrix3i Ap;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++) Ap(i,j) = std::abs(A(i,j).a % 2);
        }
        return Ap;
    }


    std::string ExactSynthesis(const U2_ZOmega& U)
    {
        auto [U_SO3, k_SO3] = U2_to_SO3(U.get_Matrix(), U.k);

        // Mat2ZOmega tmp = U.get_Matrix();
        // for(int i = 0; i < 2; i++){
        //     for(int j = 0; j < 2; j++) std::cout << tmp(i,j) << " ";
        //     std::cout << std::endl;
        // }

        Mat3ZRoot2 H_SO3, S_SO3, T_SO3;
        H_SO3(0,0) = H_SO3(0,1) = H_SO3(1,0) = H_SO3(1,2) = H_SO3(2,1) = H_SO3(2,2) = 0;
        H_SO3(0,2) = H_SO3(2,0) = 1;
        H_SO3(1,1) = -1;
        S_SO3(0,0) = S_SO3(0,2) = S_SO3(1,1) = S_SO3(1,2) = S_SO3(2,0) = S_SO3(2,1) = 0;
        S_SO3(1,0) = S_SO3(2,2) = 1;
        S_SO3(0,1) = -1;
        T_SO3(0,2) = T_SO3(1,2) = T_SO3(2,0) = T_SO3(2,1) = 0;
        T_SO3(0,0) = T_SO3(1,0) = T_SO3(1,1) = 1;
        T_SO3(0,1) = -1;
        T_SO3(2,2) = {0,1};

        Mat3ZRoot2 Hinv_SO3, Sinv_SO3, Tinv_SO3;
        Hinv_SO3 = H_SO3;
        Sinv_SO3 = S_SO3.transpose();
        Tinv_SO3 = T_SO3.transpose();

        std::string ret = "";
        while(k_SO3 > 0){
            Eigen::Matrix3i P = parity(U_SO3);

            if(P(2,0) == 0 && P(2,1) == 0 && P(2,2) == 0){
                ret += "T";
                U_SO3 = Tinv_SO3 * U_SO3;
            }else if(P(0,0) == 0 && P(0,1) == 0 && P(0,2) == 0){
                ret += "HT";
                U_SO3 = Tinv_SO3 * Hinv_SO3 * U_SO3;
            }else{
                ret += "SHT";
                U_SO3 = Tinv_SO3 * Hinv_SO3 * Sinv_SO3 * U_SO3;
            }
            k_SO3++;   // U_SO3にT_invを作用させているからインクリメント
        
            reduction(U_SO3, k_SO3);
        }

        for(auto [seq, C] : Clifford){
            bool f = true;
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++) if(U_SO3(i,j).a != C[i][j]) f = false;
            }

            if(f) ret += seq;
        }
        return ret;
    }

    int get_T_count(const U2_ZOmega& U)
    {
        auto [U_SO3, k_SO3] = U2_to_SO3(U.get_Matrix(), U.k);
        return k_SO3;
    }

    quaternion to_quaternion(std::string seq)
    {
        quaternion H(0, -1, 0, -1);
        H.unitalize();
        quaternion S(omega.real(), -omega.imag(), 0, 0);
        quaternion T(sqrt_omega.real(), -sqrt_omega.imag(), 0, 0);

        quaternion ret(1,0,0,0);
        for(auto c : seq){
            if(c == 'H') ret *= H;
            if(c == 'S') ret *= S;
            if(c == 'T') ret *= T;
        }

        return ret;
    }
}