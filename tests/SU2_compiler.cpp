#include <bits/stdc++.h>

#include "src/type.hpp"
#include "src/quaternion.hpp"
#include "src/SU2_compiler.cpp"

using namespace std;
using namespace SU2_Compiler;

int main(){
    FTYPE eps = 1e-10;
    cout << "誤差を入力してください : ";
    cin >> eps;
    quaternion U;
    // set_random_unitary_seed(1234);
    U = random_unitary();

    SU2_compiler(U, eps);
}