#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <Eigen/Core>

namespace mp = boost::multiprecision;
using namespace std;

int main(){
    mp::cpp_dec_float_100 x = 2.0f;

    // 平方根を求める
    mp::cpp_dec_float_100 result = mp::sqrt(x);

    Eigen::Matrix<mp::cpp_dec_float_100, 2, 1> v;
    v << 2.0f, 1.0f;

    double y = 2.0;

    std::cout << sqrt(x) << std::endl;
    std::cout << sqrt(y) << std::endl;

    std::cout << pow(x, x) << std::endl;
    std::cout << pow(y, y) << std::endl;
    std::cout << std::setprecision(40) << 1 / (mp::cpp_dec_float_100)3.0 << std::endl;
    std::cout << std::setprecision(40) << 1 / 3.0 << std::endl; 

    mp::cpp_dec_float_100 z = 0.5;
    double w = (double)mp::cpp_dec_float_100("1.4");
    cout << setprecision(40) << z << endl;
    cout << setprecision(40) << mp::cpp_dec_float_100("1.4") << endl;
}