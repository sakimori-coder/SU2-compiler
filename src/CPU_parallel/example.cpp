#include <bits/stdc++.h>
#include "solve_u_t.cpp"

using namespace std;

int main(){
    random_device rd;
    default_random_engine eng(1234);
    uniform_real_distribution<double> distr(0.0, 2*M_PI);

    double theta1 = distr(eng);
    double theta2 = distr(eng);
    double theta3 = distr(eng);

    complex<double> u = exp(1.0i*theta1)*cos(theta2);
    complex<double> t = exp(1.0i*theta3)*sin(theta2);

    cout << "u " << u << endl;
    cout << "t " << t << endl;
    // cout << "norm " << abs(u)*abs(u) + abs(t)*abs(t) << endl;


    tuple<ZOmega<long>, ZOmega<long>, int> ans = solve_u_t<long>(u, t, 1e-9);
    complex<double> u_prime = convert(get<0>(ans), get<2>(ans));
    complex<double> t_prime = convert(get<1>(ans), get<2>(ans));

    cout << "u' " << u_prime << endl;
    cout << "t' " << t_prime << endl;
    cout << "k " << get<2>(ans) << endl;

    cout << setprecision(15) << u + u_prime << endl;
    cout << setprecision(15) << t + t_prime <<endl;
}