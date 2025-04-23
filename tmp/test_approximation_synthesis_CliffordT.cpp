
#include <iostream>
#include <string>

#include "type.hpp"
#include "SpecialUnitary2.hpp"
#include "SU2_Compiler.hpp"

using namespace std;
using namespace Eigen;
using namespace SU2_Compiler;

int main(){
    SpecialUnitary2 V = random_unitary(123);
    Real eps = 1e-7;
    cin >> eps;
    cout << "V=" << V << endl;
    
    auto solutions = solve_approximation_synthesis_CliffordT(V, eps);

    for(auto [u, t, k] : solutions){
        SpecialUnitary2 U(u.to_Complex().real(),
                            u.to_Complex().imag(),
                            t.to_Complex().real(),
                            t.to_Complex().imag());
        U.unitalize();
        cout << U << endl;
        cout << "誤差 : " << distance(U, V) << endl;
    }

}