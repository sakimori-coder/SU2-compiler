#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <Eigen/Core>

namespace mp = boost::multiprecision;
using namespace std;

using FTYPE = mp::cpp_dec_float_50;

int main(){
    FTYPE two = 2.0;
    FTYPE sqrt2 = sqrt(two);
    cout << setprecision(50) << sqrt2 << endl;
    cout << setprecision(50) << sqrt((FTYPE)2.0) << endl;



    Eigen::Matrix<FTYPE, 2, 1> v;
    v << 1.0, 1.0;
    cout << setprecision(50) << v.norm() << endl;
    
}