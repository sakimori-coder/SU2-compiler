# test/CMakeLists.txt
file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS *.cpp)
add_executable(test_runner ${TEST_SOURCES})
add_dependencies(test_runner qd gmp mpfr mpreal eigen)

target_include_directories(test_runner PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_runner PRIVATE
    su2compiler
    QD::QD
    GMP::GMP
    GMPXX::GMPXX
    MPFR::MPFR
    MPREAL::MPREAL
    Eigen3::Eigen
    GTest::gtest_main
    Threads::Threads
)


include(GoogleTest)  # CMake >= 3.10 必須
gtest_discover_tests(test_runner
    DISCOVERY_TIMEOUT 30     # 秒
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)