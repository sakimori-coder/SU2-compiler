#include <iostream>
#include <vector>
#include <chrono>

#include "src/type.hpp"
#include "src/quaternion.hpp"
#include "src/SU2_compiler.hpp"


using namespace std;
using namespace SU2_Compiler;



int mian(){
    set_random_unitary_seed(1234);

    int n = 100;   // n個の平均を求める
    
    // ターゲットユニタリを生成
    set_random_unitary_seed(1234);
    vector<quaternion> U;
    for(int i = 0; i < n; i++){
        U[i] = random_unitary();
        cout << setprecision(20) << targetU << endl;
    }

    vector<FTYPE> eps_vec;
    FTYPE eps = 1e-1;
    while(eps >= 1e-12){
        eps_vec.push_back(eps);
        eps *= sqrt(0.1);
    }
    
    vector<double> AVE_vec;
    vector<int> MIN_vec;
    vector<int> MAX_vec;
    vector<double> AVE_TIME;
    chrono::system_clock::time_point start, end;
    double time;
    for(auto eps : eps_vec){
        cout << eps << endl;
        
        double AVE = 0.0;
        int MIN = INT_MAX;
        int MAX = 0;
        
        for(int i = 0; i < n; i++){
            start = chrono::system_clock::now();
            string V = SU2_Compiler(U[i], eps, 0);
            end = chrono::system_clock::now();
            time += static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
            int T_count = count(V.begin(), V.end(), 'T');
            AVE += (double)T_count;
            MIN = min(MIN, T_count);
            MAX = max(MAX, T_count);
            cout << "Tカウント " <<  T_count << endl;
        }
        
        AVE_vec.push_back(AVE / (double)n);
        MIN_vec.push_back(MIN);
        MAX_vec.push_back(MAX);
        AVE_TIME.push_back(time / (double)n);
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

    cout << "TIME = [";
    for(auto TIME : AVE_TIME) cout << TIME << ", ";
    cout << "]" << endl;
}