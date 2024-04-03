#include "src/Prob_Synthesis.hpp"

#include <iostream>
#include <Eigen/Core>

#include "src/type.hpp"
#include "src/quaternion.hpp"
#include "src/U2_ZOmega.hpp"
#include "src/ExactSynthesis.hpp"
#include "src/enum_u_t.hpp"

using namespace std;
using namespace SU2_Compiler;


void test_choi_matrix(){
    Eigen::Matrix<CTYPE, 2, 2> U;
    U << CTYPE(0.169466, 0.790781), CTYPE(0.527699, 0.259769),
         CTYPE(-0.477093,0.343989), CTYPE(0.035462, -0.807958);

    Eigen::Matrix<FTYPE, 4, 4> CJ_U = CJ_MB(U);

    // 多分あってる
    cout << U << endl;
    cout << CJ_U << endl;
}


void test_optimal_prob_unitary(int T_count){
    cout << "Tゲート数 : " << T_count << endl;
    quaternion targetU = random_unitary(1234);
    auto [ans, garbage] = optimal_prob_unitary(T_count, targetU, (FTYPE)0);
    FTYPE sum = 0;
    for(auto [p, U_ZOmega] : ans){
        string str = ExactSynthesis(U_ZOmega);
        cout << p << " " << str << endl; 
        sum += p;
    }
    cout << "sum " << sum << endl;
}


void test_Prob_Unitary_Synthesis(){
    FTYPE eps = 0.0001;
    quaternion targetU = random_unitary(1234);
    auto ans = Prob_Unitary_Synthesis(targetU, eps);
    for(auto [p, seq] : ans){
        cout << p << " " << seq << endl;
    }
}


void test_Random_Unitary_Synthesis(){
    const int n = 14;
    vector<quaternion> availableU(n);
    for(int i = 0; i < n; i++){
        availableU[i] = random_unitary();
    }
    
    quaternion targetU = random_unitary();

    set_targetCJUMB();
    auto [value, prob] = get_optimal_prob(availableU, targetU);
    cout << " " << value << endl;

    
    FTYPE min_v = 1.0;
    for(int i = 0; i < (1<<n); i++){
        vector<quaternion> availableU_sub;
        for(int j = 0; j < n; j++){
            if(i & (1<<j)) availableU_sub.push_back(availableU[j]);
        }
        if(availableU_sub.size() == 9){
            auto [value, prob] = get_optimal_prob(availableU_sub, targetU);
            cout << "value " << value << endl;
            min_v = min(min_v, value);
        }
    }
    cout << "min_v " << min_v << endl;
    
    // auto [value, prob] = get_optimal_prob(availableU, targetU);

    sort(prob.begin(), prob.end());
    reverse(prob.begin(), prob.end());

    FTYPE sum = 0;
    for(FTYPE p : prob){
        cout << p << endl;
        sum += p;
    }

    FTYPE min_dist = 1.0;
    for(quaternion U : availableU) min_dist = min(min_dist, distance(U, targetU));

    cout << sum << endl;
    cout << value << endl;
    cout << min_dist*min_dist << endl;
    
}

// void test_optimal_prob_unitary2(){
//     FTYPE eps = 5.29731e-09;
//     quaternion<FTYPE> targetU = random_unitary<FTYPE>(1234);
//     vector<quaternion<FTYPE>> availableU;
//     ifstream ifs("availableU.txt");
//     string line;
//     string conma;
//     while(getline(ifs, line)){
//         istringstream i_stream(line);
//         getline(i_stream, conma, ',');
//         FTYPE a(conma);
//         getline(i_stream, conma, ',');
//         FTYPE b(conma);
//         getline(i_stream, conma, ',');
//         FTYPE c(conma);
//         getline(i_stream, conma, ',');
//         FTYPE d(conma);
//         quaternion<FTYPE> U = {a,b,c,d};
//         availableU.push_back(U);
//     }

//     if(check_eps_net(availableU, targetU, eps)){
//         vector<FTYPE> prob = get_optimal_prob(availableU, targetU);
//         // cout << setprecision(30) << diamond_distance(targetU, availableU, prob) << endl;
//     }
    
// }

int main(int argc, char *argv[]){
    // test_choi_matrix();
    // test_get_optimal_prob();
    
    int T_count = 30;
    if(argc == 2) T_count = stoi(string(argv[1]));
    // test_optimal_prob_unitary(T_count);
    // test_Prob_Unitary_Synthesis();
    test_Random_Unitary_Synthesis();
}
