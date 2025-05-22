#pragma once

#include <vector>
#include "rings.hpp"

namespace su2compiler
{
    ZRoot2 pow_lambda(int n);
    std::vector<ZRoot2> one_dim_grid_problem(FTYPE x0, FTYPE x1, FTYPE y0, FTYPE y1);
}