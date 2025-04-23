#include <bits/stdc++.h>

#include "src/type.hpp"
#include "src/quaternion.hpp"
#include "src/ExactSynthesis.hpp"
#include "src/SU2_compiler.hpp"

using namespace std;
using namespace SU2_Compiler;

int main(){
    FTYPE eps = 1e-3;
    cout << "誤差を入力してください : ";
    cin >> eps;
    cout << "eps = " << eps << endl;
    quaternion U;
    set_random_unitary_seed(1234);
    U = random_unitary();
    // U = {3.0/5.0, 0, 4.0/5.0, 0.0};

    string V_str = SU2_compiler(U, eps, 0);
    cout << "V = " << V_str << endl;
    cout << "Tカウント = " << count(V_str.begin(), V_str.end(), 'T') << endl;
    quaternion V = to_quaternion(V_str);
    cout << "誤差 : " << distance(U, V) << endl;
    cout << setprecision(15) << V << endl;
}