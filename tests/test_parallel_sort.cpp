#include <bits/stdc++.h>
#include <execution>

using namespace std;

int main(){
    const int n = 7000*7000;
    vector<int> v(n);
    for(int i = 0; i < n; i++){
        v[i] = i*sin(i);
    }
    
    sort(execution::par, v.begin(), v.end());
}