#include <bits/stdc++.h>

using namespace std;


vector<int>* fun(){
    vector<int> v;
    for(int i = 0; i < 10; i++) v[i] = i;
    return &v;
}


int main(){
    vector<int>* v_prt = fun();
    vector<int> v = *v_prt;

    for(int i = 0; i < 10; i++){
        cout << v[i] << endl;
    }
}