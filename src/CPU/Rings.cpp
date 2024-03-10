#pragma once

#include <bits/stdc++.h>
#include <Eigen/Core>
// #include <boost/multiprecision/cpp_dec_float.hpp>
// #include <boost/multiprecision/cpp_int.hpp>

// using namespace std;
// namespace mp = boost::multiprecision;
// using Real_64 = mp::number<mp::cpp_dec_float<64>>;   // 仮数が64bitのfloat型

// static const Real_64 two_64 = 2.0; 
// static const Real_64 sqrt2_64 = mp::sqrt(two_64);

std::complex<double> a;

template<typename T>
struct ZRoot2
{
    T a;
    T b;

    ZRoot2() : a(0), b(0) {}
    ZRoot2(T a, T b) : a(a), b(b) {}
    ZRoot2(T a) : a(a), b(0) {}

    inline bool operator<(const ZRoot2<T> &other) const
    { 
        return a < other.a || ((a == other.a) && b < other.b);
    }

    inline bool operator>(const ZRoot2<T> &other) const
    {
        return a > other.a || ((a == other.a) && b > other.b);
    }

    inline bool operator==(const ZRoot2<T> &other) const
    {
        return a == other.a && b == other.b;
    }

    inline bool operator!=(const ZRoot2<T> &other) const
    {
        return a != other.a || b != other.b;
    }

    inline ZRoot2<T> operator+(const ZRoot2<T> &other) const
    {
        return {a + other.a, b + other.b};
    }

    inline ZRoot2<T> operator-(const ZRoot2<T> &other) const
    {
        return {a - other.a, b - other.b};
    }

    inline ZRoot2<T> operator-() const
    {
        return {-a,-b};
    }

    inline ZRoot2<T> operator*(const ZRoot2<T> &other) const
    {
        return {a*other.a + 2*b*other.b, a*other.b + b*other.a};
    }
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const ZRoot2<T>& x){
    os << "(" << x.a << ", " << x.b << ")";
    return os;
}


template<typename ITYPE, typename FTYPE>
FTYPE convert(ZRoot2<ITYPE> x){
    FTYPE sqrt2 = sqrt((FTYPE)2.0);
    return (FTYPE)x.a + sqrt2*(FTYPE)x.b;
}


// long convert_integer(ZRoot2<long> x){
//     Real_64 x_float = (Real_64)x.a + sqrt2_64*(Real_64)x.b;
//     x_float *= mp::pow(2, 60-mp::log2(x_float));
//     std::cout << (Real_64)x.a << std::endl;
//     std::cout << x_float << std::endl;
//     return (long)x_float;
// }

template<typename T>
ZRoot2<T> conj(ZRoot2<T> x){
    x.a = -x.a;
    return x;
}



/*
aω^3 + bω^2 + cω + dとして数値表現する構造体 (ω=exp(iπ/4))
*/
template<typename T>
struct ZOmega
{
    T a;
    T b;
    T c;
    T d;

    ZOmega(T A, T B, T C, T D) : a(A), b(B), c(C), d(D) {}
    ZOmega(T D) : a(0), b(0), c(0), d(D) {}
    ZOmega() : a(0), b(0), c(0), d(0) {}

    inline bool operator==(const ZOmega<T> &other) const
    {
        return (a == other.a && b == other.b) && (c == other.c && d == other.d);
    }

    inline ZOmega<T> operator+(const ZOmega<T> &other) const
    {
        return {a+other.a, b+other.b, c+other.c, d+other.d};
    }

    inline ZOmega<T> operator-(const ZOmega<T> &other) const
    {
        return {a-other.a, b-other.b, c-other.c, d-other.d};
    }

    inline ZOmega<T> operator-() const
    {
        return {-a,-b,-c,-d};
    }

    inline ZOmega<T> operator*(const ZOmega<T> &other) const
    {
        return {a*other.d + b*other.c + c*other.b + d*other.a, 
               -a*other.a + b*other.d + c*other.c + d*other.b,
               -a*other.b - b*other.a + c*other.d + d*other.c,
               -a*other.c - b*other.b - c*other.a + d*other.d};
    }

    ZOmega<T> conj(){
        return {-c, -b, -a, d};
    }

    // aω^3+bω^2+cω+d = α + iβ (α,β \in Z[√2]を仮定 <=> c+aとc-aが2で割り切れる)
    ZRoot2<T> real(){
        return {d, (c-a) / 2};
    } 
    ZRoot2<T> imag(){
        return {b, (c+a) / 2};
    }

};


template<typename T>
std::ostream& operator<<(std::ostream& os, const ZOmega<T>& x){
    os << "(" << x.a << ", " << x.b << ", " << x.c << ", " << x.d << ")";
    return os;
}


template<typename T>
ZOmega<T> std::conj(ZOmega<T> c){
    return c.conj();
}


template<typename T>
std::complex<double> convert(ZOmega<T> x){
    const std::complex<double> omega = {1.0/std::sqrt(2.0), 1.0/std::sqrt(2)};
    double a = (double)x.a;
    double b = (double)x.b;
    double c = (double)x.c;
    double d = (double)x.d;

    return a*omega*omega*omega + b*omega*omega + c*omega + d;
}


template<typename T>
std::complex<double> convert(ZOmega<T> x, int k){
    const std::complex<double> omega = {1.0/std::sqrt(2.0), 1.0/std::sqrt(2)};
    double a = (double)x.a;
    double b = (double)x.b;
    double c = (double)x.c;
    double d = (double)x.d;

    return (a*omega*omega*omega + b*omega*omega + c*omega + d) / std::pow(2.0, (double)k/2);
}


// Eigenについての動作定義


template<typename T, int N, int M>
std::ostream& operator<<(std::ostream& os, const Eigen::Matrix<ZRoot2<T>, N, M>& x){
    for(int i = 0; i < N-1; i++){
        for(int j = 0; j < M; j++) os << x(i,j) << " ";
        os << std::endl;
    }
    for(int j = 0; j < M; j++) os << x(N-1, j) << " ";
    return os;
}


template<typename T, int N, int M>
std::ostream& operator<<(std::ostream& os, const Eigen::Matrix<ZOmega<T>, N, M>& x){
    for(int i = 0; i < N-1; i++){
        for(int j = 0; j < M; j++) os << x(i,j) << " ";
        os << std::endl;
    }
    for(int j = 0; j < M; j++) os << x(N-1, j) << " ";
    return os;
}


template <typename T, int N, int M>
Eigen::Matrix<ZOmega<T>, M, N> adjoint(Eigen::Matrix<ZOmega<T>, N, M> A){
    Eigen::Matrix<ZOmega<T>, M, N> A_dag; 
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++) A_dag(j,i) = A(i,j).conj();
    }

    return A_dag;
}


template <typename T, int N, int M>
void reduction(Eigen::Matrix<ZRoot2<T>, N, M>& A, int& k){
    while(true){
        for(int i = 0; i < N; i++){
            for(int j = 0; j < M; j++){
                T a = A(i,j).a;
                T b = A(i,j).b;
                if(a % 2) return;
            }
        }

        k--;
        for(int i = 0; i < N; i++){
            for(int j = 0; j < M; j++){
            T a = A(i,j).a;
            T b = A(i,j).b;

            A(i,j).a = b;
            A(i,j).b = a / 2;
            }
        }
    }
}


template <typename T, int N, int M>
void reduction(Eigen::Matrix<ZOmega<T>, N, M>& A, int& k){
    while(true){
        for(int i = 0; i < N; i++){
            for(int j = 0; j < M; j++){
                T a = A(i,j).a;
                T b = A(i,j).b;
                T c = A(i,j).c;
                T d = A(i,j).d;
                if(a-c % 2 || b-d % 2) return;
            }
        }

        k--;
        for(int i = 0; i < N; i++){
            for(int j = 0; j < M; j++){
            T a = A(i,j).a;
            T b = A(i,j).b;
            T c = A(i,j).c;
            T d = A(i,j).d;

            A(i,j).a = (b - d) / 2;
            A(i,j).b = (c + a) / 2;
            A(i,j).c = (b + d) / 2;
            A(i,j).d = (c - a) / 2; 
            }
        }
    }
}