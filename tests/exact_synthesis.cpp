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
    string seq = Exact_synthesis(U, k);

    Eigen::Matrix<complex<double>, 2, 2> U1;
    U1 << convert(u), convert(-t.conj()), 
          convert(t), convert( u.conj());
    U1 /= pow(sqrt(2.0), k);
    cout << U1 << endl;
    
    Eigen::Matrix<complex<double>, 2, 2> H, S, T;
    H << 1., 1.,
         1.,-1.;
    H /= sqrt(2.0);
    S(0,0) = 1.;
    S(1,1) = {0.,1.};
    T(0,0) = 1.;
    T(1,1) = {1/sqrt(2.0), 1/sqrt(2.0)};

    Eigen::Matrix<complex<double>, 2, 2> U2;
    U2 << 1. , 0.,
          0. , 1.;
    for(auto c : seq){
        if(c == 'H') U2 = U2 * H;
        if(c == 'S') U2 = U2 * S;
        if(c == 'T') U2 = U2 * T;
    }

    complex<double> global_phase = U1(0,0) / U2(0,0);
    U2 = global_phase * U2;

    cout << U1 - U2 << endl;
}

bool operator<(const Eigen::Matrix<ZRoot2<int>, 3, 3>& a, const Eigen::Matrix<ZRoot2<int>, 3, 3>& b){   
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            if(a(i,j) < b(i,j)) return true;
            if(a(i,j) > b(i,j)) return false;
        }
    }
    return false;
}

void generate_Clifford(){
    Eigen::Matrix<ZOmega<int>, 2, 2> H, S;
    H(0,0) = H(0,1) = H(1,0) = {0,0,0,1};
    H(1,1) = {0,0,0,-1};
    auto [H_SO3, k_H] = to_SO3(H, 1);

    S(0,1) = S(1,0) = {0,0,0,0};
    S(0,0) = {0,0,0,1};
    S(1,1) = {0,1,0,0};
    auto [S_SO3, k_S] = to_SO3(S, 0);

    set<Eigen::Matrix<ZRoot2<int>, 3, 3>> memo;
    Eigen::Matrix<ZRoot2<int>, 3, 3> I;
    I = H_SO3 * H_SO3;
    memo.insert(I);

    int cnt = 1;

    for(int n = 1; n < 10; n++){
        for(int i = 0; i < (1<<n); i++){
            Eigen::Matrix<ZRoot2<int>, 3, 3> G = I;
            string G_str = "";
            for(int j = 0; j < n; j++){
                if(i & (1<<j)) {G = G*H_SO3; G_str += "H";}
                else {G = G*S_SO3; G_str += "S";}
            }
            
            if(memo.count(G) == 0){
                cout << "Clifford[" << cnt << "].first = \"" << G_str << "\";" << endl;
                cout << "Clifford[" << cnt << "].second << ";
                for(int x = 0; x < 3; x++){
                    for(int y = 0; y < 3; y++){
                        if(x == 2 && y == 2) continue;
                        if(x != 0 && y == 0) cout << "                  ";
                        cout << G(x,y).a << ", ";
                        if(y == 2) cout << endl;
                    }
                }
                cout << G(2,2).a << ";" << endl;

                memo.insert(G);
                cnt++;
            }
        }
    }
}

int main(){
    // test_to_SO3();
    test_Exact_synthesis();
    // generate_Clifford();

}