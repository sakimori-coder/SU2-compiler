#include <gtest/gtest.h>

#include <string>

#include "type.hpp"
#include "clifford_t/exact_synthesize.hpp"
#include "ring/Zzeta8.hpp"

using namespace std;
using namespace su2_compiler;
using su2_compiler::ring::Zzeta8;
using su2_compiler::clifford_t::U2Dzeta8;

TEST(CliffordTExactsynthesize, Ross_paper) {
    Zzeta8 u(40727366, 10614512, 10541729, -26687414);
    Zzeta8 t(20133911, 2332111, -23432014, 30805761);
    int k = 52;
    U2Dzeta8 U(u, t, 0, k);

    std::string seq1 = "HTSHTSHTSHTHTHTHTSHTHTSHTSHTSHTHTHTSHTSHTHTHTSHTHTSHTHTHTHTHTHTHTSHTSHTSHTHTSHTHTSHTHTHTHTSHTHTHTSHTHTSHTHTHTHTSHTSHTSHTHTHTSHTSHTSHTSHTHTSHTSHTSHTSHTHTSHTHTSHTSHTHTHTHTHTSHTHTHTHTSHTSHTSHTHTSHTSHTHTHTSHTHTHTHTHTSHTSHTHTHTHTHTSHTHTHTHTSHTHTHTHTHTHTH";
    std::string seq2 = clifford_t::exact_synthesize(U);
    EXPECT_EQ(seq1, seq2);
}

TEST(CliffordTExactsynthesize, Identity) {
    Zzeta8 u(1);
    Zzeta8 t(0);
    int k = 0;
    U2Dzeta8 U(u, t, 0, k);

    std::string seq1 = "";
    std::string seq2 = clifford_t::exact_synthesize(U);
    EXPECT_EQ(seq1, seq2);
}

TEST(CliffordTExactsynthesize, T) {
    Zzeta8 u(1);
    Zzeta8 t(0);
    int k = 0;
    U2Dzeta8 U(u, t, 1, k);

    std::string seq1 = "T";
    std::string seq2 = clifford_t::exact_synthesize(U);
    EXPECT_EQ(seq1, seq2);
}