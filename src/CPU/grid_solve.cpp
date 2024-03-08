#pragma once

#include <bits/stdc++.h>
#include "Rings.cpp"

template<typename T>
ZRoot2<T> get_pow_lambda(int n){
    ZRoot2<T> ans = {1,0};
    ZRoot2<T> x = {1,1};  
    if(n < 0){
        n = -n;
        x = {-1,1};
    }

    while(n > 0){
        if(n & 1) ans = ans * x;
        x = x * x;
        n >>= 1;
    }
    return ans;
}

template<typename ITYPE, typename FTYPE>
std::vector<ZRoot2<ITYPE>> one_dim_grid_problem(FTYPE x0, FTYPE x1, FTYPE y0, FTYPE y1){

    const FTYPE lambda = 1 + sqrt((FTYPE)2.0);

    ZRoot2<ITYPE> norm_factor = {1,0};
    if(x1-x0 > 1 || y1-y0 > 1){
        FTYPE n = -ceil(log(y1-y0) / (lambda));
        norm_factor = get_pow_lambda<ITYPE>(int(n));

        if(int(n)%2 == 0){
            x0 = pow(lambda,-n) * x0;
            x1 = pow(lambda,-n) * x1;
            y0 = pow(lambda, n) * y0;
            y1 = pow(lambda, n) * y1;
        }else{
            x0 = pow(lambda,-n) * x0;
            x1 = pow(lambda,-n) * x1;
            FTYPE tmp = y0;
            y0 =-pow(lambda, n) * y1;
            y1 =-pow(lambda, n) * tmp;
        }
    }

    FTYPE sqrt2 = sqrt((FTYPE)2.0);
    FTYPE sqrt8 = (FTYPE)2.0 * sqrt2;
    std::vector<ZRoot2<ITYPE>> solutions;

    // for(T a = (T)(ceil((x0+y0) / 2)); a < (T)(floor((x1+y1) / 2)) + 1; a++){
    //     for(T b = (T)(ceil((a-y1) / sqrt(2.0))); b < (T)(floor((a-y0) / sqrt(2.0))) + 1; b++){
    for(ITYPE b = (ITYPE)ceil((x0-y1) / sqrt8); b <= (ITYPE)(floor((x1-y0)) / sqrt8); b++){
        for(ITYPE a = (ITYPE)ceil(x0-(FTYPE)b*sqrt2); a <= (ITYPE)floor(x1-(FTYPE)b*sqrt2); a++){
            ZRoot2<ITYPE> candidate = {a,b};
            FTYPE cand_FTYPE = convert<ITYPE, FTYPE>(candidate);
            FTYPE conj_candi_FTYPE = convert<ITYPE, FTYPE>(conj(candidate));
            if(x0 <= cand_FTYPE && cand_FTYPE <= x1){
                if(y0 <= conj_candi_FTYPE && conj_candi_FTYPE <= y1){
                    ZRoot2 solution = norm_factor * candidate;
                    solutions.push_back(solution);
                }
            }
        }
    }

    return solutions;
} 