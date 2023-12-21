#include <bits/stdc++.h>
#include <Eigen/Dense>
#include "solve_u_t.cpp"

using namespace std;
using Eigen::Matrix2cd;


double distance(
    std::complex<double> u, std::complex<double> t, std::complex<double> u_approx, std::complex<double> t_approx)
{
    std::complex<double> diff_u = u - u_approx;
    std::complex<double> diff_t = t - t_approx;
    double ret1 =  std::sqrt((diff_u * std::conj(diff_u)).real() + (diff_t * std::conj(diff_t)).real());

    diff_u = u + u_approx;
    diff_t = t + t_approx;
    double ret2 =  std::sqrt((diff_u * std::conj(diff_u)).real() + (diff_t * std::conj(diff_t)).real());

    return std::min(ret1, ret2);
}



int main(){
    random_device rd;
    default_random_engine eng(1234);
    uniform_real_distribution<double> distr(0.0, 2*M_PI);

    vector<double> epsilons;
    double eps = sqrt(0.1);
    while(eps >= 1e-1){
        epsilons.push_back(eps);
        eps = eps / 10.0;
    }

    cout << "精度, MAX, MIN, 平均, 標準偏差" << endl;

    complex<double> exp_pi_8 = {cos(M_PI / 8.), sin(M_PI / 8.)};
    Matrix2cd T;
    T << conj(exp_pi_8), 0.,
         0., exp_pi_8;
    vector<Matrix2cd> T_powers(8);
    T_powers[0] = Matrix2cd::Identity();
    for(int i = 1; i < 8; i++) T_powers[i] = T * T_powers[i-1]; 

    for(double eps : epsilons){
        vector<int> T_counts;
        for(int i = 0; i < 100; i++){
            double theta1 = distr(eng);
            double theta2 = distr(eng);
            double theta3 = distr(eng);
            // cout << i << endl;
            complex<double> u = exp(1.0i*theta1)*cos(theta2);
            complex<double> t = exp(1.0i*theta3)*sin(theta2);

            int j = 0;
            for(j = 0; j < 8; j++){
                Matrix2cd V = T_powers[j];
                // cout << V(0,0) << " " << V(1,0) << endl;
                // cout << "distance " << distance(u, t, V(0,0), V(1,0)) << endl;
                if(distance(u, t, V(0,0), V(1,0)) < eps){
                    T_counts.push_back(j % 2);
                    cout << "aaa" << endl;
                    break;
                }
            }
            if(j != 8) continue;

            tuple<ZOmega<long>, ZOmega<long>, int> ans = solve_u_t<long>(u, t, eps);
            T_counts.push_back(2.0*(double)get<2>(ans));
            cout << 2.0*(double)get<2>(ans) << endl;
        }

        int MAX = *max_element(T_counts.begin(), T_counts.end());
        int MIN = *min_element(T_counts.begin(), T_counts.end());
        double ave = (double)accumulate(T_counts.begin(), T_counts.end(), 0) / 100.0;
        double var = 0.0;
        for(int x : T_counts) var += ((double)x - ave)*((double)x - ave) / 100.0;
        double std_dev = sqrt(var);
        
        cout << eps << ", " << MAX << ", " << MIN << ", " << ave << ", " << std_dev << endl; 
    }
}