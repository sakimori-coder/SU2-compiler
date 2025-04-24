#include <bits/stdc++.h>

#include "ring/all.hpp"

using namespace std;
using namespace su2_compiler;
using namespace ring;

int main() {
    vector<Zzeta8j> X;
    int N = 2;
    for(int a1 = -N; a1 <= N; a1++){
    for(int b1 = -N; b1 <= N; b1++){
    for(int c1 = -N; c1 <= N; c1++){
    for(int d1 = -N; d1 <= N; d1++){
        for(int a2 = -N; a2 <= N; a2++){
        for(int b2 = -N; b2 <= N; b2++){
        for(int c2 = -N; c2 <= N; c2++){
        for(int d2 = -N; d2 <= N; d2++){
            Zzeta8j U(a1,b1,c1,d1,
                      a2,b2,c2,d2);
            if(U.norm_quaternion() == Zroot2(2,0)){
                X.push_back(U);
            }
        }   
        }   
        }   
        }
    }   
    }   
    }   
    }
    cout << "|X|=" << X.size() << endl;

    vector<bool> seen(X.size(), false);
    int cnt = 0;
    for(int i = 0; i < X.size(); i++){
        if(seen[i]) continue;
        cnt++;
        cout << X[i] << endl;
        for(int j = i+1; j < X.size(); j++){
            if(X[j].leftDivisible(X[i])) seen[j] = true;
        }
    }
    cout << cnt << endl;
}