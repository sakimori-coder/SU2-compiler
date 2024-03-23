#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "src/CPU/eps_net_verification.cpp"
#include "src/CPU/quaternion.cpp"

using namespace std;
namespace mp = boost::multiprecision;


int main(){
    quaternion<mp::cpp_dec_float_50> U1(-0.542095004785001, 0.837571926662458, -0.0227331514811054, -0.0639490209293692);
    quaternion<mp::cpp_dec_float_50> U2(-0.542094660999201, 0.837571788020502, -0.0227367513122459, -0.0639524711646507);

    mp::cpp_dec_float_50 eps = 1e-5 / 2.0;
    cout << distance(U1, U2) << endl;

    cout << check_inter_2ball(U1, U2, eps) << endl;
    cout << check_inter_2ball(U2, U1, eps) << endl;
}