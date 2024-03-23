#include <iostream>

#include "src/type.hpp"
#include "src/quaternion.hpp"

using namespace std;
using namespace SU2_Compiler;

int main(){
    quaternion U = random_unitary(1234);

    cout << U << endl;
}