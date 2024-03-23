#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "src/CPU/enumerate_u_t.cpp"
#include "src/CPU/eps_net_verification.cpp"

using namespace std;
namespace mp = boost::multiprecision;
using float50 = mp::cpp_dec_float_50;
// using float50 = double;

int main(){
    random_device rd;
    default_random_engine eng(1234);
    uniform_real_distribution<double> distr_pi(0.0, M_PI);
    uniform_real_distribution<double> distr_2pi(0.0, 2*M_PI);

    double psi = distr_pi(eng);
    double theta = distr_pi(eng);
    double phi = distr_2pi(eng);
    
    double x0 = cos(psi);
    double x1 = sin(psi) * cos(theta);
    double x2 = sin(psi) * sin(theta) * cos(phi);
    double x3 = sin(psi) * sin(theta) * sin(phi);

    quaternion<float50> U(x0, x1, x2, x3);
    for(int i = 0; i < 1; i++){
        float50 eps = 1e-3;
        // float50 eps = 1e-3 * (double)i / 10.0;
        int k = 15;

        vector<quaternion<float50>> X = enumerate_u_t<long long, float50>(U, 2*eps, k);

        // for(int i = 0; i < 106; i++) X.pop_back();

        cout << "|X| = " << X.size() << endl;

        for(auto x : X) cout << distance(x, U) << endl;
        
        if(check_eps_net<float50>(X, U, eps)) cout << "OK" << endl;
        else cout << "NG" << endl;
    }
}