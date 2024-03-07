#include <bits/stdc++.h>
#include <Eigen/Core>
#include "Rings.cpp"
#include "quaternion.hpp"

// Bloch sphere representation
Eigen::Matrix<ZRoot2<T>, 3, 3> to_SO3(Eigen::Matrix<ZOmega<T>, 2, 2> U_SU2, int k_SU2){
    Eigen::Matrix<ZOmega<T>, 2, 2> U_dag;
    Eigen::Matrix<ZRoot2<T>, 3, 3> U_SO3;
    int k_SO3 = 2 * k_SU2;
    
    Eigen::Matrix<ZOmega<T>, 2, 2> X, Y, Z;
    X << {0,0,0,0}, {0,0,0,1}, 
         {0,0,0,1}, {0,0,0,0};
    Y << {0, 0,0,0}, {0,-1,0,0},
         {0,-1,0,0}, {0, 0,0,0};
    Z << {0,0,0,1}, {0,0,0,0},
         {0,0,0,0}, {0,0,0,-1};
    
    

}


template <typename T>
std::string Exact_synthesis(Eigen::Matrix<T, 2, 2> U, int k){
    
}