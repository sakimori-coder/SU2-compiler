#pragma once

#include <bits/stdc++.h>
#include "Rings.cpp"

template<typename T>
ZRoot2<T> get_pow_lambda(int n){
    ZRoot2<T> ans = {1,0};
    ZRoot2<T> x = {1,1};  
    if(n < 0){
        n = -n;
        x = {1,-1};
    }

    while(n > 0){
        if(n & 1) ans = ans * x;
        x = x * x;
        n >>= 1;
    }
    return ans;
}

template<typename T>
std::vector<ZRoot2<T>> one_dim_grid_problem(double x0, double x1, double y0, double y1){

    const double lambda = 1 + std::sqrt(2.0);

    ZRoot2<T> norm_factor = {1,0};
    if(x1-x0 > 1 || y1-y0 > 1){
        double n = -std::ceil(std::log(y1-y0) / std::log(lambda));
        norm_factor = get_pow_lambda<T>(int(n)); 

        if(int(n)%2 == 0){
            x0 = std::pow(lambda,-n) * x0;
            x1 = std::pow(lambda,-n) * x1;
            y0 = std::pow(lambda, n) * y0;
            y1 = std::pow(lambda, n) * y1;
        }else{
            x0 = std::pow(lambda,-n) * x0;
            x1 = std::pow(lambda,-n) * x1;
            double tmp = y0;
            y0 =-std::pow(lambda, n) * y1;
            y1 =-std::pow(lambda, n) * tmp;
        }
    }

    std::vector<ZRoot2<T>> solutions;
    
    for(T a = (T)(std::ceil((x0+y0) / 2)); a < (T)(std::floor((x1+y1) / 2)) + 1; a++){
        for(T b = (T)(std::ceil((a-y1) / std::sqrt(2.0))); b < (T)(std::floor((a-y0) / std::sqrt(2.0))) + 1; b++){
            ZRoot2<T> candidate = {a,b};
            double cand_double = convert(candidate);
            double conj_candi_double = convert(conj(candidate));
            if(x0 <= cand_double && cand_double <= x1){
                if(y0 <= conj_candi_double && conj_candi_double <= y1){
                    ZRoot2 solution = norm_factor * candidate;
                    solutions.push_back(solution);
                }
            }
        }
    }

    return solutions;
} 