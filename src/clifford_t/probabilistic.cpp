#include "clifford_t/probabilistic.hpp"

#include "type.hpp"
#include "probabilistic_tools/mixing_su2.hpp"
#include "clifford_t/exact_synthesis.hpp"
#include "clifford_t/deterministic.hpp"

namespace su2compiler::clifford_t
{


std::vector<std::pair<Real, std::string>> probabilistic_synthesis(SU2 V, Real eps){
    if(eps >= Real(1.0)) {
        return {{1.0, ""}};
    }

    int t = 0;
    while(true){
        auto res = probabilistic_synthesis_fixed_t(V, eps, t);
        if(!res.empty()) return res;
        t++;
    }

    return {};   // Failure
}

    
std::vector<std::pair<Real, std::string>> probabilistic_synthesis_fixed_t(SU2 V, Real eps, UINT t){
    Real delta = pow(Real(2.0), -t/3);

    while(true) {
        std::vector<U2Dzeta8> availableU;
        for(int s = 0; s <= t; s++) {
            auto availableU_s = deterministic_synthesis_fixed_t(V, min(2*delta, Real(1.0)), s);
            availableU.insert(availableU.end(), availableU_s.begin(), availableU_s.end());
        }
        
        std::vector<SU2> availableU_SU2;
        for(auto U : availableU) availableU_SU2.push_back(U.to_SU2());
        
        if(2*delta >= Real(1.0) || availableU.size() > 50) {
            auto [opt_distance, opt_distribution] = mixing_su2::optimize_distribution(availableU_SU2, V);
            
            if(opt_distance <= eps) {
                Real tol = 1e-15;
                std::vector<std::pair<Real, std::string>> ret;
                for(int i = 0; i < availableU.size(); i++){
                    if(opt_distribution[i] < tol) continue;
                    ret.push_back({opt_distribution[i], exact_synthesis(availableU[i])});
                }
                return ret;
            }

            delta += delta / 10;
        }
    }
}



}

