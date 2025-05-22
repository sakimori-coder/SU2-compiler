cmake_minimum_required(VERSION 3.24)
project(su2compiler LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(LIB_DIR ${CMAKE_SOURCE_DIR}/lib)

set(USE_TEST Yes)

include(ExternalProject)
include(FetchContent)

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)



#----------------------------------------------------------------------
# 外部モジュールのbuild
#----------------------------------------------------------------------
set(EXTERN_PREFIX ${CMAKE_BINARY_DIR}/extern)
file(MAKE_DIRECTORY ${EXTERN_PREFIX}/include ${EXTERN_PREFIX}/lib)

##### GMP 6.3.0 #####
ExternalProject_Add(
    gmp
    URL "https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz"
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    PREFIX ${CMAKE_BINARY_DIR}/gmp
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ./configure --enable-static
                                  --enable-pthreads
                                  --disable-shared 
                                  --with-pic
                                  --enable-cxx
                                  --prefix=${EXTERN_PREFIX}
    BUILD_COMMAND     make -j
    INSTALL_COMMAND   make install
    BUILD_BYPRODUCTS ${EXTERN_PREFIX}/lib/libgmp.a ${EXTERN_PREFIX}/lib/libgmpxx.a
)
# GMP
add_library(GMP::GMP STATIC IMPORTED)
add_dependencies(GMP::GMP gmp)
set_target_properties(GMP::GMP PROPERTIES
    IMPORTED_LOCATION             ${EXTERN_PREFIX}/lib/libgmp.a
    INTERFACE_INCLUDE_DIRECTORIES ${EXTERN_PREFIX}/include)

# GMPXX
add_library(GMPXX::GMPXX STATIC IMPORTED)
add_dependencies(GMPXX::GMPXX gmp)
set_target_properties(GMPXX::GMPXX PROPERTIES
    IMPORTED_LOCATION             ${EXTERN_PREFIX}/lib/libgmpxx.a
    INTERFACE_LINK_LIBRARIES      GMP::GMP
    INTERFACE_INCLUDE_DIRECTORIES ${EXTERN_PREFIX}/include)


##### MPFR 4.2.2 #####
ExternalProject_Add(
    mpfr
    URL "https://www.mpfr.org/mpfr-current/mpfr-4.2.2.tar.xz"
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    PREFIX ${CMAKE_BINARY_DIR}/mpfr
    DEPENDS gmp
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ./configure --enable-static
                                  --enable-thread-safe
                                  --disable-shared
                                  --with-pic
                                  --with-gmp=${EXTERN_PREFIX}
                                  --prefix=${EXTERN_PREFIX}
    BUILD_COMMAND     make -j
    INSTALL_COMMAND   make install
    BUILD_BYPRODUCTS ${EXTERN_PREFIX}/lib/libmpfr.a
)
add_library(MPFR::MPFR STATIC IMPORTED)
add_dependencies(MPFR::MPFR mpfr)
set_target_properties(MPFR::MPFR PROPERTIES
    IMPORTED_LOCATION             ${EXTERN_PREFIX}/lib/libmpfr.a
    INTERFACE_INCLUDE_DIRECTORIES ${EXTERN_PREFIX}/include
    INTERFACE_LINK_LIBRARIES      GMP::GMP
)


##### MPREAL 3.7.0 #####
ExternalProject_Add(
    mpreal
    GIT_REPOSITORY  https://github.com/advanpix/mpreal.git
    GIT_TAG         mpfrc++-3.7.0
    PREFIX          ${CMAKE_BINARY_DIR}/mpreal
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND 
        ${CMAKE_COMMAND} -E copy
        <SOURCE_DIR>/mpreal.h
        ${EXTERN_PREFIX}/include/mpreal.h          
)
add_library(MPREAL::MPREAL STATIC IMPORTED)
add_dependencies(MPREAL::MPREAL mpreal)
set_target_properties(MPREAL::MPREAL PROPERTIES
    IMPORTED_LOCATION             ${EXTERN_PREFIX}/lib/libmpfr.a
    INTERFACE_INCLUDE_DIRECTORIES ${EXTERN_PREFIX}/include
    INTERFACE_LINK_LIBRARIES      MPFR::MPFR
)


##### EIGEN 3.4.0 #####
ExternalProject_Add(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    PREFIX ${CMAKE_BINARY_DIR}/eigen
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy_directory 
            <SOURCE_DIR>/Eigen 
            ${EXTERN_PREFIX}/include/Eigen
    TEST_COMMAND ""
)
add_library(Eigen3::Eigen INTERFACE IMPORTED)
add_dependencies(Eigen3::Eigen eigen)
set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${EXTERN_PREFIX}/include
    SYSTEM TRUE)


##### GTEST 1.12.1 #####
if(USE_TEST)
    enable_testing()
    
    find_package(GTest QUIET)

    if(NOT TARGET GTest::gtest)
        set(FETCHCONTENT_TRY_FIND_PACKAGE OFF)
        set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

        FetchContent_Declare(
        googletest_fetch
        GIT_REPOSITORY https://github.com/google/googletest
        GIT_TAG release-1.12.1
        )
        FetchContent_MakeAvailable(googletest_fetch) 
    endif()
else()
	message(STATUS "Skip downloding googletest")
endif()



#----------------------------------------------------------------------
# SRCコードのbuild
#----------------------------------------------------------------------
add_subdirectory(src)


#----------------------------------------------------------------------
# TESTコードのbuild
#----------------------------------------------------------------------
if(USE_TEST)
    add_subdirectory(test)
endif()

