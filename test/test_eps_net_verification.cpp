#include <bits/stdc++.h>
#include <Eigen/Core>
#include "../src/CPU/eps_net_verification.cpp"

using namespace std;


void test_check_inter_2ball(){
    cout << "start test_check_inter2ball" << endl;
    random_device rd;
    default_random_engine eng(rd()); 
    uniform_real_distribution<double> distr_pi(0.0, M_PI);
    uniform_real_distribution<double> distr_2pi(0.0, 2*M_PI);

    double psi = distr_pi(eng);
    double theta = distr_pi(eng);
    double phi = distr_2pi(eng);
    
    double x0 = cos(psi);
    double x1 = sin(psi) * cos(theta);
    double x2 = sin(psi) * sin(theta) * cos(phi);
    double x3 = sin(psi) * sin(theta) * sin(phi);

    psi -= 0.001;
    double y0 = cos(psi);
    double y1 = sin(psi) * cos(theta);
    double y2 = sin(psi) * sin(theta) * cos(phi);
    double y3 = sin(psi) * sin(theta) * sin(phi);

    quaternion<double> u1(x0, x1, x2, x3);
    quaternion<double> u2(y0, y1, y2, y3);

    double d = distance(u1, u2);   // u1とu2の直線距離
    double eps = d / 2.0;   // この大きさのε-ball同士は接する

    // 交わらないはず
    if(check_inter_2ball(u1, u2, eps - 1e-5)) cout << "交わる" << endl;
    else cout << "交わらない" << endl;

    // 交わるはず
    if(check_inter_2ball(u1, u2, eps + 1e-5)) cout << "交わる" << endl;
    else cout << "交わらない" << endl;

}


void test_cal_inter_3ball(){
    cout << "start test_call_inter3ball" << endl;
    random_device rd;
    default_random_engine eng(4202536160);   // 交点を持つ乱数に設定した
    uniform_real_distribution<double> distr_pi(0.0, M_PI);
    uniform_real_distribution<double> distr_2pi(0.0, 2*M_PI);

    double psi = distr_pi(eng);
    double theta = distr_pi(eng);
    double phi = distr_2pi(eng);
    
    double x0 = cos(psi);
    double x1 = sin(psi) * cos(theta);
    double x2 = sin(psi) * sin(theta) * cos(phi);
    double x3 = sin(psi) * sin(theta) * sin(phi);

    psi -= 0.00001;
    double y0 = cos(psi);
    double y1 = sin(psi) * cos(theta);
    double y2 = sin(psi) * sin(theta) * cos(phi);
    double y3 = sin(psi) * sin(theta) * sin(phi);

    theta -= 0.00001;
    double z0 = cos(psi);
    double z1 = sin(psi) * cos(theta);
    double z2 = sin(psi) * sin(theta) * cos(phi);
    double z3 = sin(psi) * sin(theta) * sin(phi);

    quaternion<double> u1(x0, x1, x2, x3);
    quaternion<double> u2(y0, y1, y2, y3);
    quaternion<double> u3(z0, z1, z2, z3);

    double eps = 0.01;
    auto [inter1, inter2] = cal_inter_3ball(u1, u2, u3, eps);

    // 全てεになっているはず
    cout << "d(u1,inter1) = " << distance(u1, inter1) << endl;   
    cout << "d(u2,inter1) = " << distance(u2, inter1) << endl;
    cout << "d(u3,inter1) = " << distance(u3, inter1) << endl;
    cout << "d(u1,inter2) = " << distance(u1, inter2) << endl;   
    cout << "d(u2,inter2) = " << distance(u2, inter2) << endl;
    cout << "d(u3,inter2) = " << distance(u3, inter2) << endl;
    // 交点のノルムは1になっているはず
    cout << "||inter1|| = " << inter1.norm() << endl;
    cout << "||inter2|| = " << inter2.norm() << endl;

    double d = distance(u1, u2);
}

int main(){
    test_cal_inter_3ball();
    test_check_inter_2ball();
}