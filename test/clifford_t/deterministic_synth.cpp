#include <gtest/gtest.h>

#include "core/type.hpp"
#include "core/su2.hpp"
#include "ring/all.hpp"
#include "clifford_t/deterministic_synth.hpp"

using namespace su2compiler;


template <typename Real>
class CliffordTDetSynthTest : public ::testing::Test {
protected:
    void SetUp() override {
        if constexpr(std::is_same_v<Real, double>) {
            exit(0);
        }
        if constexpr(std::is_same_v<Real, mpfr::mpreal>) {
            mpfr::mpreal::set_default_prec(256);
        }
    }
};
using RealTypes = ::testing::Types<REAL_SCALAR_TYPE_LIST>;
TYPED_TEST_SUITE(CliffordTDetSynthTest, RealTypes);


TYPED_TEST(CliffordTDetSynthTest, random) {
    using RealType = TypeParam;

    SU2<RealType> V = random_unitary<RealType>(1234);
    RealType eps = 1e-5;

    std::string sequence = clifford_t::deterministic::synth(V, eps);
    // std::cout << "T-count = " << std::count(sequence.begin(), sequence.end(), 'T') << std::endl;

    SU2<RealType> H(Complex<RealType>(0, 1) / math::SQRT2<RealType>(), 
                    Complex<RealType>(0, 1) / math::SQRT2<RealType>());
    SU2<RealType> S(std::conj(math::ZETA8<RealType>()), Complex<RealType>(0, 0));
    SU2<RealType> T(std::conj(math::ZETA16<RealType>()), Complex<RealType>(0, 0));

    SU2<RealType> U(1,0,0,0);
    for(char gate : sequence) {
        switch (gate) {
            case 'H':
                U *= H;
                break;
            case 'S':
                U *= S;
                break;
            case 'T':
                U *= T;
                break;
            default:
                break;
        }
    }

    RealType dist = distance(U, V);
    EXPECT_TRUE(dist < eps);
}