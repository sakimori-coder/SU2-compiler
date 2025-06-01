#include <gtest/gtest.h>

#include <string>

#include "core/type.hpp"
#include "ring/Zzeta8.hpp"
#include "clifford_t/exact_synth.hpp"

using namespace std;
using namespace su2compiler;
using su2compiler::ring::Zzeta8;
using su2compiler::clifford_t::exact::U2Dzeta8;

TEST(CliffordT_ExactSynth_Test, Ross_paper) {
    Zzeta8 u(40727366, 10614512, 10541729, -26687414);
    Zzeta8 t(20133911, 2332111, -23432014, 30805761);
    Natural k = 52;
    U2Dzeta8 U(u, t, 0, k);

    std::string seq1 = "HTSHTSHTSHTHTHTHTSHTHTSHTSHTSHTHTHTSHTSHTHTHTSHTHTSHTHTHTHTHTHTHTSHTSHTSHTHTSHTHTSHTHTHTHTSHTHTHTSHTHTSHTHTHTHTSHTSHTSHTHTHTSHTSHTSHTSHTHTSHTSHTSHTSHTHTSHTHTSHTSHTHTHTHTHTSHTHTHTHTSHTSHTSHTHTSHTSHTHTHTSHTHTHTHTHTSHTSHTHTHTHTHTSHTHTHTHTSHTHTHTHTHTHTH";
    std::string seq2 = clifford_t::exact::synth(U);
    EXPECT_EQ(seq1, seq2);
}

TEST(CliffordT_ExactSynth_Test, Identity) {
    Zzeta8 u(1);
    Zzeta8 t(0);
    Natural k = 0;
    U2Dzeta8 U(u, t, 0, k);

    std::string seq1 = "";
    std::string seq2 = clifford_t::exact::synth(U);
    EXPECT_EQ(seq1, seq2);
}

TEST(CliffordT_ExactSynth_Test, T) {
    Zzeta8 u(1);
    Zzeta8 t(0);
    Natural k = 0;
    U2Dzeta8 U(u, t, 1, k);

    std::string seq1 = "T";
    std::string seq2 = clifford_t::exact::synth(U);
    EXPECT_EQ(seq1, seq2);
}