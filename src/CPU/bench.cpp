#include <bits/stdc++.h>
#include "solve_u_t.cpp"

using namespace std;

int main(){
    random_device rd;
    default_random_engine eng(1234);
    uniform_real_distribution<double> distr(0.0, 2*M_PI);

    vector<double> epsilons;
    double eps = sqrt(0.1);
    while(eps >= 1e-9){
        epsilons.push_back(eps);
        eps = eps / 10.0;
    }

    cout << "精度, MAX, MIN, 平均, 標準偏差" << endl;

    for(double eps : epsilons){
        vector<int> T_counts;
        for(int i = 0; i < 100; i++){
            double theta1 = distr(eng);
            double theta2 = distr(eng);
            double theta3 = distr(eng);
            // cout << i << endl;
            complex<double> u = exp(1.0i*theta1)*cos(theta2);
            complex<double> t = exp(1.0i*theta3)*sin(theta2);

            tuple<ZOmega<long>, ZOmega<long>, int> ans = solve_u_t<long>(u, t, eps);
            T_counts.push_back(2.0*(double)get<2>(ans));
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