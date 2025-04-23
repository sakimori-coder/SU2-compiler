#pragma once

#include <string>
#include "type.hpp"
#include "quaternion.hpp"



namespace SU2_Compiler{
    std::string SU2_compiler(quaternion U, FTYPE eps, int l);
}