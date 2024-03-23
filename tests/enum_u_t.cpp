#include <bits/stdc++.h>
#include "src/type.hpp"
#include "src/quaternion.hpp"
#include "src/U2_ZOmega.hpp"
#include "src/enum_u_t.hpp"

using namespace std;
using namespace SU2_Compiler;

int main(){
    FTYPE eps = 0.0001;
    quaternion targetU = random_unitary(1234);
    cout << targetU << endl;

    auto availableU = enum_u_t(targetU, eps, 20, {});

    cout << availableU.first.size() << endl;
}