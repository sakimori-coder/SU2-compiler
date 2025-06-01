#include <gtest/gtest.h>

#include "core/type.hpp"
#include "math/number_theory.hpp"

using namespace su2compiler;

bool isPrimeBF(const Integer& n){
    if(n == 1) return false;
    for(Integer i = 2; i*i <= n; i++){
        if(n % i == 0) return false;
    }
    return true;
}


TEST(NumberTheoryTest, isPrime)
{
    using math::number_theory::isPrime;

    for(Integer n = 1000; n < 1100; n++) {
        EXPECT_TRUE(isPrime(n) == isPrimeBF(n));
    }

    Integer n50("47373862195627799992848703730035240754888050152689");
    EXPECT_TRUE(isPrime(n50));

}