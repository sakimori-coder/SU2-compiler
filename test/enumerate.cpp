#include <bits/stdc++.h>
#include "src/CPU/quaternion.hpp"
#include "src/CPU/enumerate_u_t.cpp"

using namespace std;

using FTYPE = double;

int main(){
    FTYPE eps = 0.0001;
    quaternion<FTYPE> targetU = random_unitary<FTYPE>(1234);

    vector<quaternion<FTYPE>> availableU = enumerate_u_t<long long, FTYPE>(targetU, eps, 20);

    cout << availableU.size() << endl;
}