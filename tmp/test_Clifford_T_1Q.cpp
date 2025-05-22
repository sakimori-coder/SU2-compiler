#include <cassert>
#include <iostream>
#include <string>

#include "type.hpp"
#include "rings.hpp"
#include "SpecialUnitary2.hpp"
#include "Clifford_T_1Q.hpp"

using namespace std;
using namespace Eigen;
using namespace su2compiler;


int main(){
    ZOmega u(40727366, 10614512, 10541729, -26687414);
    ZOmega t(20133911, 2332111, -23432014, 30805761);
    Integer k = 52;

    Clifford_T_1Q U(u, t, 0, k);
    string seq1 = U.compute_sequence();
    
    string seq2 = "HTSHTSHTSHTHTHTHTSHTHTSHTSHTSHTHTHTSHTSHTHTHTSHTHTSHTHTHTHTHTHTHTSHTSHTSHTHTSHTHTSHTHTHTHTSHTHTHTSHTHTSHTHTHTHTSHTSHTSHTHTHTSHTSHTSHTSHTHTSHTSHTSHTSHTHTSHTHTSHTSHTHTHTHTHTSHTHTHTHTSHTSHTSHTHTSHTSHTHTHTSHTHTHTHTHTSHTSHTHTHTHTHTSHTHTHTHTSHTHTHTHTHTHTH";

    assert(seq1 == seq2);

    Real theta = PI / 128;
    Real eps = 1e-10;
    SpecialUnitary2 V(cos(-theta/2), sin(-theta/2), 0.0, 0.0);
    
    assert(distance(U.to_SU2(), V) < eps);
}