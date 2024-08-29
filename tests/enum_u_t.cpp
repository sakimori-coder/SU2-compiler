#include <bits/stdc++.h>
#include "src/type.hpp"
#include "src/quaternion.hpp"
#include "src/U2_ZOmega.hpp"
#include "src/enum_u_t.hpp"
#include "src/ExactSynthesis.hpp"

using namespace std;
using namespace SU2_Compiler;

int main(){
    FTYPE eps = pow((FTYPE)2.0, -(FTYPE)16.0);
    int k = 15;
    FTYPE pre_eps = 0.0;
    set_random_unitary_seed(1234);
    // quaternion targetU = random_unitary();

    FTYPE three("3.0");
    quaternion targetU(sqrt(2.0), 1.0 + sqrt(2.0), 1.0 - sqrt(2.0), three);
    targetU.unitalize();
    
    // cout << setprecision(15) << targetU << "\n" << endl;

    std::vector<ZRoot2> XY;
    std::vector<ZRoot2> ZW;
    std::vector<U2_ZOmega> candidates;

    for(int k = 1; k < 40; k++){
        cout << "k = " << k << endl;
        auto availableU = enum_u_t(targetU, eps, pre_eps, k, 0, XY, ZW, candidates);

        cout << "|X| = " << availableU.size() << endl;

        // for(auto U_omega : availableU){
        //     cout << "Tカウント " << get_T_count(U_omega) << endl;
        //     cout << to_quaternion(U_omega) << endl;
        // }
    
    }
}