file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS *.cpp)

add_library(su2compiler STATIC ${SRC_FILES})

set_target_properties(su2compiler PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

target_include_directories(su2compiler PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_compile_options(su2compiler PUBLIC -O3 -fdiagnostics-show-template-tree -fno-elide-type)
target_link_libraries(su2compiler PUBLIC
    GMP::GMP
    GMPXX::GMPXX
    MPFR::MPFR
    MPREAL::MPREAL
    Eigen3::Eigen
    Threads::Threads
    OpenMP::OpenMP_CXX
)