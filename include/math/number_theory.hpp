#pragma once

#include <random>

#include "core/type.hpp"


namespace su2compiler::math::number_theory
{

Integer pow_mod(
    Integer x, 
    Integer y,
    const Integer& mod
)
{
    Integer ret(1);
    x %= mod;
    while(y > 0) {
        if(y % 2 == 1) {
            ret = (ret * x) % mod; 
        }
        x = (x * x) % mod;
        y >>= 1;
    }
    return ret;
}

bool isPrime(const Integer& n, int k = 20)
{
    if(n < 2) return false;
    if(n == 2 || n == 3) return true;
    if(n % 2 == 0) return false;
    
    Integer d = n-1;
    int s = 0;
    while(d % 2 == 0) {
        d >>= 1;
        s++;
    }

    std::random_device rd;
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(rd());
    for(int i = 0; i < k; i++) {
        Integer a = rng.get_z_range(n-4) + 2;
        Integer x = pow_mod(a, d, n);
        if(x == 1 || x == n-1) continue;

        bool witness = true;
        for(int r = 1; r < s; r++) {
            x = (x * x) % n;
            if(x == n-1) {
                witness = false;
                break;
            }
        }

        if(witness) return false;
    }
    return true;
}

}