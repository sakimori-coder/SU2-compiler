#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace mp = boost::multiprecision;

template<typename T>
T sqrt(T x) {return mp::sqrt(x);}

double sqrt(double x) {return std::sqrt(x);}