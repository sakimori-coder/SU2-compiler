#include <bits/stdc++.h>
#include <Eigen/Core>
#include "src/CPU/quaternion.hpp"

using namespace std;

int main(){
    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<double> distr_pi(0.0, M_PI);
    uniform_real_distribution<double> distr_2pi(0.0, 2*M_PI);

    double psi = distr_pi(eng);
    double theta = distr_pi(eng);
    double phi = distr_2pi(eng);
    
    double x0 = cos(psi);
    double x1 = sin(psi) * cos(theta);
    double x2 = sin(psi) * sin(theta) * cos(phi);
    double x3 = sin(psi) * sin(theta) * sin(phi);

    psi -= 0.001;
    double y0 = cos(psi);
    double y1 = sin(psi) * cos(theta);
    double y2 = sin(psi) * sin(theta) * cos(phi);
    double y3 = sin(psi) * sin(theta) * sin(phi);

    quaternion<double> u(x0, x1, x2, x3);
    quaternion<double> v(y0, y1, y2, y3);

    cout << distance(u, v) << endl;

    Eigen::Matrix2cd U = u.get_Matrix();
    Eigen::Matrix2cd V = v.get_Matrix();
    cout << U << endl;
    cout << V << endl;
}
