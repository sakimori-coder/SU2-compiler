#include "clifford_t/exact_synth.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "core/type.hpp"
#include "core/su2.hpp"
#include "ring/all.hpp"


namespace su2compiler::clifford_t::exact {

static const U2Dzeta8 H_U2(ring::Zzeta8(1), ring::Zzeta8(1), 4, 1);
static const U2Dzeta8 S_U2(ring::Zzeta8(1), ring::Zzeta8(0), 2, 0);
static const U2Dzeta8 Sinv_U2 = S_U2.adjoint();
static const U2Dzeta8 T_U2(ring::Zzeta8(1), ring::Zzeta8(0), 1, 0);
static const U2Dzeta8 Tinv_U2 = T_U2.adjoint();

static const SO3Droot2 H_SO3 = H_U2.toSO3Droot2();
static const SO3Droot2 S_SO3 = S_U2.toSO3Droot2();
static const SO3Droot2 Sinv_SO3 = S_SO3.transpose();
static const SO3Droot2 T_SO3 = T_U2.toSO3Droot2();
static const SO3Droot2 Tinv_SO3 = T_SO3.transpose();


//==============================================================================
//  implement SO3Droot2
//==============================================================================
SO3Droot2::SO3Droot2(Matrix3Zroot2 _mat, Natural _k) {
	mat = _mat;
	k = _k;
	// if(mat * mat.transpose() != ((Integer)1 << k) * Matrix3Zroot2::Identity()){
	// 	throw std::invalid_argument("SO3Droot2() : the given matrix is not in SO(3)");
	// }
	this->reduce();
}

void SO3Droot2::reduce() {
	ring::Zroot2 sqrt2(0,1);
	while(k >= 0) {
		bool divisible = true;
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) if(!mat(i,j).divisible(sqrt2)) divisible = false;
		}

		if(!divisible) break;
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) mat(i,j) /= sqrt2;
		}
		k--;
	}
}

SO3Droot2& SO3Droot2::operator*=(const SO3Droot2& r) {
	mat *= r.mat;
	k += r.k;
	this->reduce();
	return *this;
}



//==============================================================================
//  implement U2Zzeta8
//==============================================================================
U2Dzeta8::U2Dzeta8(Matrix2Zzeta8 _mat, Natural _k) {
	mat = _mat;
	k = _k;
	Matrix2Zzeta8 mat_dag;
	mat_dag << mat(0,0).conj_complex(), mat(1,0).conj_complex(),
			   mat(0,1).conj_complex(), mat(1,1).conj_complex();
	if(mat * mat_dag != ring::Zzeta8((Integer(1) << _k)) * Matrix2Zzeta8::Identity()) {
		throw std::invalid_argument("U2Dzeta8() : the given matrix is not in U(2)");
	}
	this->reduce();
}

U2Dzeta8::U2Dzeta8(ring::Zzeta8 u, ring::Zzeta8 t, int l, Natural _k) {
	if(u.norm_complex() + t.norm_complex() != ring::Zroot2(Integer(1) << _k)) {
		throw std::invalid_argument("U2Dzeta8() : the given matrix is not in U(2)");
	}

	l %= 8;
	ring::Zzeta8 phase(1);
	for(int i = 0; i < l; i++) phase *= ring::Zzeta8(0,1,0,0);
	mat << u, -t.conj_complex() * phase,
		   t,  u.conj_complex() * phase;
	k = _k;
	this->reduce();
}

U2Dzeta8::U2Dzeta8(std::string sequence) {
	mat = Matrix2Zzeta8::Identity();
	k = 0;

	for(char gate : sequence) {
		switch (gate)
		{
		case 'H':
			*this *= H_U2;
			break;
		case 'S':
			*this *= S_U2;
			break;
		case 'T':
			*this *= T_U2;
			break;
		default:
			throw std::invalid_argument(std::string("U2Dzeta8() : Unsupported gate: '") + gate + "'");
		}
	}
}

template <typename RealType>
[[nodiscard]] SU2<RealType> U2Dzeta8::toSU2() const {
	int l;
	ring::Zzeta8 zeta8_pow(1);
	for(l = 0; l < 8; l++) {
		if(mat(1,1) / mat(0,0).conj_complex() == zeta8_pow) break;
		zeta8_pow *= ring::Zzeta8(0,1,0,0); 
	}

	std::complex<RealType> u = mat(0,0).toComplex<RealType>() / math::pow_ui(math::SQRT2<RealType>(), k);
	std::complex<RealType> t = mat(1,0).toComplex<RealType>() / math::pow_ui(math::SQRT2<RealType>(), k);
	for(int i = 0; i < l; i++) {
		u /= math::ZETA16<RealType>();
		t /= math::ZETA16<RealType>();
	}
	return SU2<RealType>(u, t);
}
#define X(RealType) template [[nodiscard]] SU2<RealType> U2Dzeta8::toSU2() const;
REAL_SCALAR_TYPES
#undef X


[[nodiscard]] SO3Droot2 U2Dzeta8::toSO3Droot2() const {
	Matrix2Zzeta8 mat_dag;
	mat_dag << mat(0,0).conj_complex(), mat(1,0).conj_complex(),
			   mat(0,1).conj_complex(), mat(1,1).conj_complex();

	Matrix2Zzeta8 X,Y,Z;
	X << ring::Zzeta8(0), ring::Zzeta8(1),
		 ring::Zzeta8(1), ring::Zzeta8(0);
	Y << ring::Zzeta8(0), ring::Zzeta8(0,0,-1,0),
		 ring::Zzeta8(0,0,1,0), ring::Zzeta8(0);
	Z << ring::Zzeta8(1), ring::Zzeta8(0),
		 ring::Zzeta8(0), ring::Zzeta8(-1);

    Matrix2Zzeta8 UXUdag, UYUdag, UZUdag;
    UXUdag = mat * X * mat_dag * ring::Zzeta8(2);
    UYUdag = mat * Y * mat_dag * ring::Zzeta8(2);
    UZUdag = mat * Z * mat_dag * ring::Zzeta8(2);

	Matrix3Zroot2 mat_SO3;
	int k_SO3 = 2*k + 2;  // The “+2” is there because we multiply U(Pauli)U† by 2.

	auto trace_real = [](Matrix2Zzeta8 M) { 
		ring::Zzeta8 trace = M(0,0) + M(1,1);
		return ring::Zroot2(trace.a / 2, (trace.b - trace.d) / 4);	
	};
    
	mat_SO3(0,0) = trace_real(X * UXUdag);
    mat_SO3(0,1) = trace_real(X * UYUdag);
    mat_SO3(0,2) = trace_real(X * UZUdag);

    mat_SO3(1,0) = trace_real(Y * UXUdag);
    mat_SO3(1,1) = trace_real(Y * UYUdag);
    mat_SO3(1,2) = trace_real(Y * UZUdag);

    mat_SO3(2,0) = trace_real(Z * UXUdag);
    mat_SO3(2,1) = trace_real(Z * UYUdag);
    mat_SO3(2,2) = trace_real(Z * UZUdag);

	return SO3Droot2(mat_SO3, k_SO3);
}

void U2Dzeta8::reduce() {
	ring::Zroot2 sqrt2(0,1);
	while(k >= 0) {
		bool divisible = true;
		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 2; j++) {
				if(!mat(i,j).divisible(sqrt2)) divisible = false;
			}
		}

		if(!divisible) break;
		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 2; j++) mat(i,j) /= sqrt2;
		}
		k--;
	}
}

U2Dzeta8& U2Dzeta8::operator*=(const U2Dzeta8& r) {
	mat *= r.mat;
	k += r.k;
	this->reduce();
	return *this;
}




extern const std::pair<std::string, Matrix3Zroot2> Clifford_SO3_list[24];

std::string synth(U2Dzeta8 U) {
	return synth(U.toSO3Droot2());
}

std::string synth(SO3Droot2 U) {
    std::string sequence = "";
    while(U.k > 0){
        Eigen::Matrix3i P;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++) P(i,j) = ((U.mat(i,j).a & 1) != 0);
        }


        if(P(2,0) == 0 && P(2,1) == 0 && P(2,2) == 0){
            sequence += "T";
			U = Tinv_SO3 * U;
        }else if(P(0,0) == 0 && P(0,1) == 0 && P(0,2) == 0){
            sequence += "HT";
            U = Tinv_SO3 * H_SO3 * U;
        }else{
            sequence += "SHT";
            U = Tinv_SO3 * H_SO3 * Sinv_SO3 * U;
        }
    }

    for(auto [Clifford_seq, Clifford_SO3] : Clifford_SO3_list){
        if((U.mat.array() == Clifford_SO3.array()).all()) {
			sequence += Clifford_seq;
		}
	}

	return sequence;
}

Natural get_Tcount(const U2Dzeta8& U) {
	SO3Droot2 U_SO3 = U.toSO3Droot2();
	return U_SO3.k;
}


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