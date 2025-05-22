#include "Clifford_T_1Q.hpp"

#include <string>
#include <utility>
#include <complex>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Core>

#include "type.hpp"
#include "rings.hpp"
#include "SpecialUnitary2.hpp"

namespace su2compiler
{

template<typename Matrix>
static void reduce(Matrix& M, Integer& k)
{
    int n = M.rows();
    while(k >= 0){  
        bool dividable = true;
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(!M(i,j).divisibleBySqrt2()) dividable = false; 
            }
        }
        if(!dividable) break;

        k--;
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                M(i,j).divideBySqrt2();
            }
        }
    }
}

Clifford_T_1Q::Clifford_T_1Q()
{
    mat << 1, 0,
           0, 1;
    k = 0;
}

Clifford_T_1Q::Clifford_T_1Q(ZOmega u, ZOmega t, int l, Integer k)
{     
    if(u * complex_conj(u) + t * complex_conj(t) != ((Integer)1<<k)){
        throw std::invalid_argument(std::string("Clifford_T_1Q: norm of (u,t)^t must be 1"));
    }
    l %= 8;
    mat << u, -complex_conj(t) * omega_pow[l],
           t,  complex_conj(u) * omega_pow[l];
    this->k = k;
}

Clifford_T_1Q::Clifford_T_1Q(const std::string& sequence)
{
    Matrix2ZOmega Hmat, Smat, Tmat;
    Hmat << 1, 1,
            1,-1;
    Smat << 1, 0,
            0, ZOmega(0,0,1,0);
    Tmat << 1, 0,
            0, ZOmega(0,1,0,0);
    
    mat << 1, 0,
           0, 1;
    k = 0;
    for(char gate : sequence){
        switch (gate) {
            case 'H':
                mat *= Hmat;
                k++;
                break;
            case 'S':
                mat *= Smat;
                break;
            case 'T':
                mat *= Tmat;
                break;
            default:
                throw std::invalid_argument(std::string("Unsupported gate: '") + gate + "'");
        }
    }

    
    reduce(mat, k);
}


extern const Matrix2ZOmega X_ZOmega;
extern const Matrix2ZOmega Y_ZOmega;
extern const Matrix2ZOmega Z_ZOmega;

extern const Matrix3ZRoot2 Hinv_SO3;
extern const Matrix3ZRoot2 Sinv_SO3;
extern const Matrix3ZRoot2 Tinv_SO3;

extern const std::pair<std::string, Matrix3ZRoot2> Clifford_SO3[24];


std::string Clifford_T_1Q::compute_sequence()
{
    ZOmega sqrt2(0,1,0,-1);

    Matrix2ZOmega U = this->mat;
    U *= sqrt2;
    Matrix2ZOmega U_dag;
    U_dag << complex_conj(U(0,0)), complex_conj(U(1,0)),
             complex_conj(U(0,1)), complex_conj(U(1,1));

    Matrix2ZOmega UXU_dag, UYU_dag, UZU_dag;
    UXU_dag = U * X_ZOmega * U_dag;
    UYU_dag = U * Y_ZOmega * U_dag;
    UZU_dag = U * Z_ZOmega * U_dag;

    Matrix3ZRoot2 U_SO3;
    Integer k_SO3 = 2*(this->k) + 2;   // The reason of +2 is U *= sqrt2

    U_SO3(0,0) = UXU_dag(0,1).real();
    U_SO3(0,1) = UYU_dag(0,1).real();
    U_SO3(0,2) = UZU_dag(0,1).real();

    U_SO3(1,0) = -UXU_dag(0,1).imag();
    U_SO3(1,1) = -UYU_dag(0,1).imag();
    U_SO3(1,2) = -UZU_dag(0,1).imag();

    U_SO3(2,0) = UXU_dag(0,0).real();
    U_SO3(2,1) = UYU_dag(0,0).real();
    U_SO3(2,2) = UZU_dag(0,0).real();

    reduce(U_SO3, k_SO3);

    std::string ret = "";
    while(k_SO3 > 0){
        Eigen::Matrix3i P;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++) P(i,j) = (Integer)std::abs(int(U_SO3(i,j).a % 2));
        }

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
        reduce(U_SO3, k_SO3);
    }

    for(auto [seq, CliffordMat] : Clifford_SO3){
        if((U_SO3.array() == CliffordMat.array()).all()) ret += seq;
    }

    return ret;
}


int Clifford_T_1Q::compute_Tcount() const
{
    ZOmega sqrt2(-1,0,1,0);

    Matrix2ZOmega U = this->mat;
    U *= sqrt2;
    Matrix2ZOmega U_dag;
    U_dag << complex_conj(U(0,0)), complex_conj(U(1,0)),
             complex_conj(U(0,1)), complex_conj(U(1,1));

    Matrix2ZOmega UXU_dag, UYU_dag, UZU_dag;
    UXU_dag = U * X_ZOmega * U_dag;
    UYU_dag = U * Y_ZOmega * U_dag;
    UZU_dag = U * Z_ZOmega * U_dag;

    Matrix3ZRoot2 U_SO3;
    Integer k_SO3 = 2*(this->k) + 2;   // The reason of +2 is U *= sqrt2

    U_SO3(0,0) = UXU_dag(0,1).real();
    U_SO3(0,1) = UYU_dag(0,1).real();
    U_SO3(0,2) = UZU_dag(0,1).real();

    U_SO3(1,0) = -UXU_dag(0,1).imag();
    U_SO3(1,1) = -UYU_dag(0,1).imag();
    U_SO3(1,2) = -UZU_dag(0,1).imag();

    U_SO3(2,0) = UXU_dag(0,0).real();
    U_SO3(2,1) = UYU_dag(0,0).real();
    U_SO3(2,2) = UZU_dag(0,0).real();

    reduce(U_SO3, k_SO3);

    return k_SO3;
}


SpecialUnitary2 Clifford_T_1Q::to_SU2() const
{
    Complex phase = {1.0, 0.0};
    for(int i = 0; i < 8; i++){
        if(complex_conj(mat(0,0)) * omega_pow[i] == mat(1,1)) break;
        phase *= sqrt_omega;
    }
    
    Complex u = mat(0,0).to_Complex() * std::conj(phase);
    Complex t = mat(1,0).to_Complex() * std::conj(phase);

    u /= pow(sqrt2, (Real)k);
    t /= pow(sqrt2, (Real)k);

    return SpecialUnitary2(u.real(), u.imag(), t.real(), t.imag());
}



ZOmega omega_pow[8] = {{1,0,0,0}, {0,1,0,0}, {0,0,1,0},  {0,0,0,1},
                       {-1,0,0,0}, {0,-1,0,0}, {0,0,-1,0}, {0,0,0,-1}};


const Matrix2ZOmega X_ZOmega = (Matrix2ZOmega() <<
        0, 1,
        1, 0
).finished();

const Matrix2ZOmega Y_ZOmega = (Matrix2ZOmega() <<
        0, ZOmega(0,0,-1,0),
        ZOmega(0,0,1,0), 0
).finished();

const Matrix2ZOmega Z_ZOmega = (Matrix2ZOmega() <<
        1, 0,
        0,-1
).finished();



const Matrix3ZRoot2 Hinv_SO3 = (Matrix3ZRoot2() <<
        0, 0, 1,
        0,-1, 0,
        1, 0, 0
).finished();

const Matrix3ZRoot2 Sinv_SO3 = (Matrix3ZRoot2() <<
        0, 1, 0,
       -1, 0, 0,
        0, 0, 1
).finished();

const Matrix3ZRoot2 Tinv_SO3 = (Matrix3ZRoot2() <<
        1, 1, 0,
       -1, 1, 0,
        0, 0, ZRoot2(0,1)
).finished();



const std::pair<std::string, Matrix3ZRoot2> Clifford_SO3[24] = 
{
    {"", (Matrix3ZRoot2() << 
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
    ).finished()},

    {"S", (Matrix3ZRoot2() << 
            0,-1, 0,
            1, 0, 0,
            0, 0, 1
    ).finished()},

    {"H", (Matrix3ZRoot2() << 
            0, 0, 1,
            0,-1, 0,
            1, 0, 0
    ).finished()},

    {"SS", (Matrix3ZRoot2() << 
           -1, 0, 0,
            0,-1, 0,
            0, 0, 1
    ).finished()},

    {"HS", (Matrix3ZRoot2() << 
            0, 0, 1,
           -1, 0, 0,
            0,-1, 0
    ).finished()},

    {"SH", (Matrix3ZRoot2() << 
            0, 1, 0,
            0, 0, 1,
            1, 0, 0
    ).finished()},

    {"SSS", (Matrix3ZRoot2() << 
            0, 1, 0,
           -1, 0, 0,
            0, 0, 1
    ).finished()},

    {"HSS", (Matrix3ZRoot2() << 
            0, 0, 1,
            0, 1, 0,
           -1, 0, 0
    ).finished()},

    {"SHS", (Matrix3ZRoot2() << 
            1, 0, 0,
            0, 0, 1,
            0,-1, 0
    ).finished()},

    {"SSH", (Matrix3ZRoot2() << 
            0, 0,-1,
            0, 1, 0,
            1, 0, 0
    ).finished()},

    {"HSH", (Matrix3ZRoot2() << 
            1, 0, 0,
            0, 0,-1,
            0, 1, 0
    ).finished()},

    {"HSSS", (Matrix3ZRoot2() << 
            0, 0, 1,
            1, 0, 0,
            0, 1, 0
    ).finished()},

    {"SHSS", (Matrix3ZRoot2() << 
            0,-1, 0,
            0, 0, 1,
           -1, 0, 0
    ).finished()},

    {"SSHS", (Matrix3ZRoot2() << 
            0, 0,-1,
            1, 0, 0,
            0,-1, 0
    ).finished()},

    {"HSHS", (Matrix3ZRoot2() << 
            0,-1, 0,
            0, 0,-1,
            1, 0, 0
    ).finished()},

    {"HSSH", (Matrix3ZRoot2() << 
            1, 0, 0,
            0,-1, 0,
            0, 0,-1
    ).finished()},

    {"SHSSS", (Matrix3ZRoot2() << 
           -1, 0, 0,
            0, 0, 1,
            0, 1, 0
    ).finished()},

    {"SSHSS", (Matrix3ZRoot2() << 
            0, 0,-1,
            0,-1, 0,
           -1, 0, 0
    ).finished()},

    {"HSHSS", (Matrix3ZRoot2() << 
           -1, 0, 0,
            0, 0,-1,
            0,-1, 0
    ).finished()},

    {"HSSHS", (Matrix3ZRoot2() << 
            0,-1, 0,
           -1, 0, 0,
            0, 0,-1
    ).finished()},

    {"SHSSH", (Matrix3ZRoot2() << 
            0, 1, 0,
            1, 0, 0,
            0, 0,-1
    ).finished()},

    {"SSHSSS", (Matrix3ZRoot2() << 
            0, 0,-1,
           -1, 0, 0,
            0, 1, 0
    ).finished()},

    {"HSHSSS", (Matrix3ZRoot2() << 
            0, 1, 0,
            0, 0,-1,
           -1, 0, 0
    ).finished()},

    {"HSSHSS", (Matrix3ZRoot2() << 
           -1, 0, 0,
            0, 1, 0,
            0, 0,-1
    ).finished()},
};
    
}