#pragma once

#include <bits/stdc++.h>
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

    ZRoot2() {}
    ZRoot2(T a, T b) : a(a), b(b) {}
    ZRoot2(T a) : a(a), b(0) {}

    inline bool operator<(const ZRoot2<T> &other) const
    { 
        // if(a == other.a) return (b < other.b);
        // return (a < other.a);
        return a < other.a || (!(other.a < a) && b < other.b);
    }

    inline bool operator==(const ZRoot2<T> &other) const
    {
        return a == other.a && b == other.b;
    }

    inline ZRoot2<T> operator+(const ZRoot2<T> &other) const
    {
        return {a + other.a, b + other.b};
    }

    inline ZRoot2<T> operator-(const ZRoot2<T> &other) const
    {
        return {a - other.a, b - other.b};
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


template<typename T>
double convert(ZRoot2<T> x){
    return x.a + sqrt(2.0)*x.b;
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

    inline bool operator==(const ZOmega<T> &other) const
    {
        return (a == other.a && b == other.b) && (c == other.c && d == other.d);
    }

};


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