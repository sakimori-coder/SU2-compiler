#include <bits/stdc++.h>
#include <Eigen/Core>

#include "src/CPU/Rings.cpp"
#include "src/CPU/exact_decomp.cpp"

using namespace std;

void test_to_SO3(){
    Eigen::Matrix<ZOmega<int>, 2, 2> T, H, S;
    T(0,0) = {0,0,0,1};
    T(0,1) = {0,0,0,0};
    T(1,0) = {0,0,0,0};
    T(1,1) = {0,0,1,0};
    auto [T_SO3, k_T] = to_SO3(T, 0);
    cout << k_T << endl;
    cout << T_SO3 << endl;

    H(0,0) = H(0,1) = H(1,0) = {0,0,0,1};
    H(1,1) = {0,0,0,-1};
    auto [H_SO3, k_H] = to_SO3(H, 1);
    cout << k_H << endl;
    cout << H_SO3 << endl;

    S(0,1) = S(1,0) = {0,0,0,0};
    S(0,0) = {0,0,0,1};
    S(1,1) = {0,1,0,0};
    auto [S_SO3, k_S] = to_SO3(S, 0);
    cout << k_S << endl;
    cout << S_SO3 << endl;
}

void test_Exact_synthesis(){
    ZOmega<long long> u = {190, 214, 595, -582};
    ZOmega<long long> t = {91, 305, -20, -415};
    int k = 20;

    cout << u*u.conj() + t*t.conj() << endl;
    cout << (1<<20) << endl;

    Eigen::Matrix<ZOmega<long long>, 2, 2> U;
    U << u, -t.conj(),
         t,  u.conj();
    Exact_synthesis(U, k); 
}

int main(){
    // test_to_SO3();
    test_Exact_synthesis();
    
}