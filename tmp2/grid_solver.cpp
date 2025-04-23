
#include "grid_solver.hpp"

#include <iostream>
#include <vector>
#include "type.hpp"
#include "rings.hpp"


namespace SU2_Compiler
{
    ZRoot2 pow_lambda(int n){
        ZRoot2 ans = {1,0};
        ZRoot2 x = {1,1};   // nが正ならλ, 負ならλ^-1
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

    std::vector<ZRoot2> one_dim_grid_problem(FTYPE x0, FTYPE x1, FTYPE y0, FTYPE y1){
        const FTYPE lambda = 1 + sqrt2;

        ZRoot2 norm_factor = {1,0};
        if(x1-x0 > 1 || y1-y0 > 1){
            FTYPE n = -ceil(log(y1-y0) / (lambda));
            norm_factor = pow_lambda(int(n));

            // if(std::abs(int(n)) % 2 == 0){
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
        std::vector<ZRoot2> solutions;

        FTYPE eps = 1e-20;   // 数値誤差のため
        for(ITYPE b = (ITYPE)ceil((x0-y1) / sqrt8); b <= (ITYPE)floor((x1-y0) / sqrt8); b++){
            for(ITYPE a = (ITYPE)ceil(x0-(FTYPE)b*sqrt2); a <= (ITYPE)floor(x1-(FTYPE)b*sqrt2); a++){
                ZRoot2 candidate = {a,b};
                FTYPE cand_FTYPE = ZRoot2_to_FTYPE(candidate);
                ZRoot2 conj_candidate = conj(candidate);
                FTYPE conj_candi_FTYPE = ZRoot2_to_FTYPE(conj_candidate);
                if(x0 - eps <= cand_FTYPE && cand_FTYPE <= x1 + eps){
                    if(y0 - eps <= conj_candi_FTYPE && conj_candi_FTYPE <= y1 + eps){
                        ZRoot2 solution = norm_factor * candidate;
                        solutions.push_back(solution);
                    }
                }
            }
        }

        return solutions;
    } 
}