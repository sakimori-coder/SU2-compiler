#include "clifford_t/exact_synthesis.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "type.hpp"
#include "ring/Zroot2.hpp"
#include "ring/Zzeta8.hpp"


namespace su2_compiler {
namespace clifford_t {

extern const Matrix2Zzeta8 X_Zzeta8;
extern const Matrix2Zzeta8 Y_Zzeta8;
extern const Matrix2Zzeta8 Z_Zzeta8;

extern const Matrix3Zroot2 Hinv_SO3;
extern const Matrix3Zroot2 Sinv_SO3;
extern const Matrix3Zroot2 Tinv_SO3;

extern const std::pair<std::string, Matrix3Zroot2> Clifford_SO3_list[24];


std::string exact_synthesis(ring::Zzeta8 u, ring::Zzeta8 t, int l, Integer k) {
    if(u.norm_complex() + t.norm_complex() != ring::Zroot2(Integer(1) << k, 0)){
        std::ostringstream oss;
        oss << "clifford_t::exact_synthesis : required |u|² + |t|² = " + std::to_string(Integer(1) << k)
            << ", but |u|² + |t|² = " << u.norm_complex() + t.norm_complex();
        throw std::domain_error(oss.str());
    }

    Matrix2Zzeta8 U, U_dag;
    if(l % 2 == 0){
        U << u, -t.conj_complex(),
             t,  u.conj_complex();        
    }else{
        U << u, -t.conj_complex() * ring::Zzeta8(0,1,0,0),
             t,  u.conj_complex() * ring::Zzeta8(0,1,0,0);
    }
    U *= ring::Zzeta8(ring::Zroot2(0,1));

    U_dag << U(0,0).conj_complex(), U(1,0).conj_complex(),
             U(0,1).conj_complex(), U(1,1).conj_complex();

    Matrix2Zzeta8 UXU_dag, UYU_dag, UZU_dag;
    UXU_dag = U * X_Zzeta8 * U_dag;
    UYU_dag = U * Y_Zzeta8 * U_dag;
    UZU_dag = U * Z_Zzeta8 * U_dag; 

    Matrix3Zroot2 U_SO3;
    Integer k_SO3 = 2 * k + 2;

    auto get_real = [](ring::Zzeta8& x) { return ring::Zroot2(x.a, (x.b-x.d) / 2); };
    auto get_imag = [](ring::Zzeta8& x) { return ring::Zroot2(x.c, (x.b+x.d) / 2); };
    U_SO3(0,0) = get_real(UXU_dag(0,1));
    U_SO3(0,1) = get_real(UYU_dag(0,1));
    U_SO3(0,2) = get_real(UZU_dag(0,1));

    U_SO3(1,0) = -get_imag(UXU_dag(0,1));
    U_SO3(1,1) = -get_imag(UYU_dag(0,1));
    U_SO3(1,2) = -get_imag(UZU_dag(0,1));

    U_SO3(2,0) = get_real(UXU_dag(0,0));
    U_SO3(2,1) = get_real(UYU_dag(0,0));
    U_SO3(2,2) = get_real(UZU_dag(0,0));

    auto reduce = [](Matrix3Zroot2 &U, Integer& k){
        while(k > 0){
            if(U.array().unaryExpr([](ring::Zroot2 x){ return x.divisible(ring::Zroot2(0,1)); }).all()){
                U /= ring::Zroot2(0,1);
                k--;
            }else{
                break;
            }
        }
    };
    reduce(U_SO3, k_SO3);
    
    std::string sequence = "";
    while(k_SO3 > 0){
        Eigen::Matrix3i P;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++) P(i,j) = ((U_SO3(i,j).a & 1) != 0);
        }

        if(P(2,0) == 0 && P(2,1) == 0 && P(2,2) == 0){
            sequence += "T";
            U_SO3 = Tinv_SO3 * U_SO3;
        }else if(P(0,0) == 0 && P(0,1) == 0 && P(0,2) == 0){
            sequence += "HT";
            U_SO3 = Tinv_SO3 * Hinv_SO3 * U_SO3;
        }else{
            sequence += "SHT";
            U_SO3 = Tinv_SO3 * Hinv_SO3 * Sinv_SO3 * U_SO3;
        }
        k_SO3++;   // U_SO3にT_invを作用させているからインクリメント
        reduce(U_SO3, k_SO3);
    }

    for(auto [seq, Clifford_SO3] : Clifford_SO3_list){
        if((U_SO3.array() == Clifford_SO3.array()).all()) sequence += seq;
    }

    return sequence;
}



const Matrix2Zzeta8 X_Zzeta8 = (Matrix2Zzeta8() <<
        ring::Zzeta8(0), ring::Zzeta8(1),
        ring::Zzeta8(1), ring::Zzeta8(0)
).finished();

const Matrix2Zzeta8 Y_Zzeta8 = (Matrix2Zzeta8() <<
        ring::Zzeta8(0), ring::Zzeta8(0,0,-1,0),
        ring::Zzeta8(0,0,1,0), ring::Zzeta8(0)
).finished();

const Matrix2Zzeta8 Z_Zzeta8 = (Matrix2Zzeta8() <<
        ring::Zzeta8(1), ring::Zzeta8(0),
        ring::Zzeta8(0), ring::Zzeta8(-1)
).finished();




const Matrix3Zroot2 Hinv_SO3 = (Matrix3Zroot2() <<
        0, 0, 1,
        0,-1, 0,
        1, 0, 0
).finished();

const Matrix3Zroot2 Sinv_SO3 = (Matrix3Zroot2() <<
        0, 1, 0,
       -1, 0, 0,
        0, 0, 1
).finished();

const Matrix3Zroot2 Tinv_SO3 = (Matrix3Zroot2() <<
        1, 1, 0,
       -1, 1, 0,
        0, 0, ring::Zroot2(0,1)
).finished();




const std::pair<std::string, Matrix3Zroot2> Clifford_SO3_list[24] = 
{
    {"", (Matrix3Zroot2() << 
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
    ).finished()},

    {"S", (Matrix3Zroot2() << 
            0,-1, 0,
            1, 0, 0,
            0, 0, 1
    ).finished()},

    {"H", (Matrix3Zroot2() << 
            0, 0, 1,
            0,-1, 0,
            1, 0, 0
    ).finished()},

    {"SS", (Matrix3Zroot2() << 
           -1, 0, 0,
            0,-1, 0,
            0, 0, 1
    ).finished()},

    {"HS", (Matrix3Zroot2() << 
            0, 0, 1,
           -1, 0, 0,
            0,-1, 0
    ).finished()},

    {"SH", (Matrix3Zroot2() << 
            0, 1, 0,
            0, 0, 1,
            1, 0, 0
    ).finished()},

    {"SSS", (Matrix3Zroot2() << 
            0, 1, 0,
           -1, 0, 0,
            0, 0, 1
    ).finished()},

    {"HSS", (Matrix3Zroot2() << 
            0, 0, 1,
            0, 1, 0,
           -1, 0, 0
    ).finished()},

    {"SHS", (Matrix3Zroot2() << 
            1, 0, 0,
            0, 0, 1,
            0,-1, 0
    ).finished()},

    {"SSH", (Matrix3Zroot2() << 
            0, 0,-1,
            0, 1, 0,
            1, 0, 0
    ).finished()},

    {"HSH", (Matrix3Zroot2() << 
            1, 0, 0,
            0, 0,-1,
            0, 1, 0
    ).finished()},

    {"HSSS", (Matrix3Zroot2() << 
            0, 0, 1,
            1, 0, 0,
            0, 1, 0
    ).finished()},

    {"SHSS", (Matrix3Zroot2() << 
            0,-1, 0,
            0, 0, 1,
           -1, 0, 0
    ).finished()},

    {"SSHS", (Matrix3Zroot2() << 
            0, 0,-1,
            1, 0, 0,
            0,-1, 0
    ).finished()},

    {"HSHS", (Matrix3Zroot2() << 
            0,-1, 0,
            0, 0,-1,
            1, 0, 0
    ).finished()},

    {"HSSH", (Matrix3Zroot2() << 
            1, 0, 0,
            0,-1, 0,
            0, 0,-1
    ).finished()},

    {"SHSSS", (Matrix3Zroot2() << 
           -1, 0, 0,
            0, 0, 1,
            0, 1, 0
    ).finished()},

    {"SSHSS", (Matrix3Zroot2() << 
            0, 0,-1,
            0,-1, 0,
           -1, 0, 0
    ).finished()},

    {"HSHSS", (Matrix3Zroot2() << 
           -1, 0, 0,
            0, 0,-1,
            0,-1, 0
    ).finished()},

    {"HSSHS", (Matrix3Zroot2() << 
            0,-1, 0,
           -1, 0, 0,
            0, 0,-1
    ).finished()},

    {"SHSSH", (Matrix3Zroot2() << 
            0, 1, 0,
            1, 0, 0,
            0, 0,-1
    ).finished()},

    {"SSHSSS", (Matrix3Zroot2() << 
            0, 0,-1,
           -1, 0, 0,
            0, 1, 0
    ).finished()},

    {"HSHSSS", (Matrix3Zroot2() << 
            0, 1, 0,
            0, 0,-1,
           -1, 0, 0
    ).finished()},

    {"HSSHSS", (Matrix3Zroot2() << 
           -1, 0, 0,
            0, 1, 0,
            0, 0,-1
    ).finished()},
};


}
}