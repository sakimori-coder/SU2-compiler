#include <bits/stdc++.h>
#include <Eigen/Core>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "src/CPU/Prob_synthesis.cpp"
#include "src/CPU/enumerate_u_t.cpp"
#include "src/CPU/exact_synthesis.cpp"

using namespace std;
namespace mp = boost::multiprecision;

using FTYPE = mp::cpp_dec_float_50;
// using FTYPE = double;

void test_choi_matrix(){
    Eigen::Matrix<complex<FTYPE>, 2, 2> U;
    U << complex<FTYPE>(0.169466, 0.790781), complex<FTYPE>(0.527699, 0.259769),
         complex<FTYPE>(-0.477093,0.343989), complex<FTYPE>(0.035462, -0.807958);

    Eigen::Matrix<FTYPE, 4, 4> CJ_U = CJ_MB<FTYPE>(U);

    // 多分あってる
    cout << U << endl;
    cout << CJ_U << endl;
}


// void test_get_optimal_prob(){
//     random_device rd;
//     default_random_engine eng(1234);
//     uniform_real_distribution<double> distr_pi(0.0, M_PI);
//     uniform_real_distribution<double> distr_2pi(0.0, 2*M_PI);

//     FTYPE psi = distr_pi(eng);
//     FTYPE theta = distr_pi(eng);
//     FTYPE phi = distr_2pi(eng);
    
//     FTYPE x0 = cos(psi);
//     FTYPE x1 = sin(psi) * cos(theta);
//     FTYPE x2 = sin(psi) * sin(theta) * cos(phi);
//     FTYPE x3 = sin(psi) * sin(theta) * sin(phi);

//     quaternion<FTYPE> U(x0, x1, x2, x3);
//     FTYPE eps = 0.9*1e-3;
//     cout << "eps^2 = " << eps*eps << endl;
//     int k = 18;
//     auto [X, pre_results] = enumerate_u_t<long long, FTYPE>(U, eps, k);

//     std::vector<FTYPE> prob = get_optimal_prob(X, U);

//     FTYPE sum = 0;
//     for(auto p : prob){
//         cout << p << " ";
//         sum += p;
//     }
//     cout << endl;

//     cout << diamond_distance<FTYPE>(U, X, prob) << endl;
// }


void test_optimal_prob_unitary(int T_count){
    cout << "Tゲート数 : " << T_count << endl;
    quaternion<FTYPE> targetU = random_unitary<FTYPE>(1234);
    auto [ans, garbage] = optimal_prob_unitary(T_count, targetU);
    FTYPE sum = 0;
    for(auto [p, unitary] : ans){
        auto [u, t, k] = unitary;
        string str = Exact_synthesis(u, t, k, T_count % 2);
        cout << p << " " << str << endl; 
        sum += p;
    }
    cout << "sum " << sum << endl;
}

void test_optimal_prob_unitary2(){
    FTYPE eps = 5.29731e-09;
    quaternion<FTYPE> targetU = random_unitary<FTYPE>(1234);
    vector<quaternion<FTYPE>> availableU;
    ifstream ifs("availableU.txt");
    string line;
    string conma;
    while(getline(ifs, line)){
        istringstream i_stream(line);
        getline(i_stream, conma, ',');
        FTYPE a(conma);
        getline(i_stream, conma, ',');
        FTYPE b(conma);
        getline(i_stream, conma, ',');
        FTYPE c(conma);
        getline(i_stream, conma, ',');
        FTYPE d(conma);
        quaternion<FTYPE> U = {a,b,c,d};
        availableU.push_back(U);
    }

    if(check_eps_net(availableU, targetU, eps)){
        vector<FTYPE> prob = get_optimal_prob(availableU, targetU);
        // cout << setprecision(30) << diamond_distance(targetU, availableU, prob) << endl;
    }
    
}

int main(int argc, char *argv[]){
    // test_choi_matrix();
    // test_get_optimal_prob();
    
    int T_count = 30;
    if(argc == 2) T_count = stoi(string(argv[1]));
    test_optimal_prob_unitary(T_count);

    // test_optimal_prob_unitary2();
    
}
