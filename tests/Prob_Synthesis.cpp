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
    auto [ans, garbage] = optimal_prob_unitary(T_count, targetU, 0);
    FTYPE sum = 0;
    for(auto [p, U_ZOmega] : ans){
        string str = ExactSynthesis(U_ZOmega);
        cout << p << " " << str << endl; 
        sum += p;
    }
    cout << "sum " << sum << endl;
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
    test_optimal_prob_unitary(T_count);

    // test_optimal_prob_unitary2();
    
}
