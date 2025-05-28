#include <gtest/gtest.h>

#include "type.hpp"
#include "su2.hpp"
#include "probabilistic_tools/mixing_su2.hpp"

using namespace::su2compiler;

TEST(mixing_su2, Choi_Jamiolkowski) {
    SU2 U = su2compiler::random_unitary();
    Matrix4C S = mixing_su2::Choi_Jamiolkowski(U.to_EigenMatrix());

    EXPECT_NEAR(double( (S - S.adjoint()).cwiseAbs().sum() ), 0.0, 1e-10);
    EXPECT_NEAR(double( std::abs((S*S).trace()) ), 1.0, 1e-10);
}

TEST(mixing_su2, Choi_Jamiolkowski_MagicBasis) {
    SU2 U = su2compiler::random_unitary();
    Matrix4C S = mixing_su2::Choi_Jamiolkowski_MagicBasis(U.to_EigenMatrix());

    EXPECT_NEAR(double( (S - S.transpose()).cwiseAbs().sum() ), 0.0, 1e-10);
    EXPECT_NEAR(double( std::abs((S*S).trace()) ), 1.0, 1e-10);
}


TEST(mixing_su2, optimize_distribution) {
    mpfr::mpreal::set_default_prec(300);
    SU2 V = su2compiler::random_unitary();
    int N = 50;
    std::vector<SU2> availableU(N);
    for(int i = 0; i < N; i++) availableU[i] = su2compiler::random_unitary();

    for(int i = 0; i < N; i++){
        std::cout << distance(V, availableU[i]) << std::endl;
    }

    auto [opt_distance, opt_distribution] = mixing_su2::optimize_distribution(availableU, V);

    std::cout << opt_distance << std::endl;
}