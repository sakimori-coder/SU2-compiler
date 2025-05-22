#pragma once

#include <string>
#include "type.hpp"
#include "quaternion.hpp"



namespace su2compiler{
    std::string su2compiler(quaternion U, FTYPE eps, int l);
}