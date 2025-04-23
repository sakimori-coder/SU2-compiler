#include <gtest/gtest.h>

#include <string>

#include "type.hpp"
#include "clifford_t/exact_synthesis.hpp"
#include "ring/Zzeta8.hpp"

using namespace std;
using namespace su2_compiler;
using su2_compiler::ring::Zzeta8;

TEST(ct_exact_synthesis, ross_paper) {
    Zzeta8 u(40727366, 10614512, 10541729, -26687414);
    Zzeta8 t(20133911, 2332111, -23432014, 30805761);
    Integer k = 52;

    std::string seq1 = "HTSHTSHTSHTHTHTHTSHTHTSHTSHTSHTHTHTSHTSHTHTHTSHTHTSHTHTHTHTHTHTHTSHTSHTSHTHTSHTHTSHTHTHTHTSHTHTHTSHTHTSHTHTHTHTSHTSHTSHTHTHTSHTSHTSHTSHTHTSHTSHTSHTSHTHTSHTHTSHTSHTHTHTHTHTSHTHTHTHTSHTSHTSHTHTSHTSHTHTHTSHTHTHTHTHTSHTSHTHTHTHTHTSHTHTHTHTSHTHTHTHTHTHTH";
    std::string seq2 = clifford_t::exact_synthesis(u, t, 0, k);
    EXPECT_EQ(seq1, seq2);
}