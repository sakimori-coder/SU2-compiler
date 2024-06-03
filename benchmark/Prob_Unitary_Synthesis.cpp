#include <iostream>
#include <vector>

#include "src/type.hpp"
#include "src/quaternion.hpp"
#include "src/Prob_Synthesis.hpp"

using namespace std;
using namespace SU2_Compiler;

int main(){
    set_random_unitary_seed(1234);
    
    int n = 100;   // n個の平均を求める
    
    vector<FTYPE> eps_vec;
    FTYPE eps = 0.1;
    while(eps > 1e-5){
        eps_vec.push_back(eps);
        eps *= sqrt(0.1);
    }
    
    vector<double> AVE_vec;
    vector<int> MIN_vec;
    vector<int> MAX_vec;
    for(auto eps : eps_vec){
        cout << eps << endl;
        
        double AVE = 0.0;
        int MIN = INT_MAX;
        int MAX = 0;
        
        for(int i = 0; i < n; i++){
            quaternion targetU = random_unitary();
            cout << setprecision(20) << targetU << endl;
            auto mixed_unitary = Prob_Unitary_Synthesis(targetU, eps);
            int T_count = get_T_count(mixed_unitary);
            AVE += (double)T_count;
            MIN = min(MIN, T_count);
            MAX = max(MAX, T_count);
            cout << "Tカウント " <<  T_count << endl;
        }
        
        AVE_vec.push_back(AVE / n);
        MIN_vec.push_back(MIN);
        MAX_vec.push_back(MAX);
    }
    
    cout << "eps = [";
    for(auto eps : eps_vec) cout << eps << ", ";
    cout << "]" << endl;

    cout << "AVE = [";
    for(auto AVE : AVE_vec) cout << AVE << ", ";
    cout << "]" << endl;

    cout << "MIN = [";
    for(auto MIN : MIN_vec) cout << MIN << ", ";
    cout << "]" << endl;

    cout << "MAX = [";
    for(auto MAX : MAX_vec) cout << MAX << ", ";
    cout << "]" << endl;
}