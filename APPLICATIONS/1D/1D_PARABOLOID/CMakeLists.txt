MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Project name and version
project(shallow_water)

# Enable Fortran
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 ")
    set(bounds "-Wall -fbounds-check -fimplicit-none -fcheck=all")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE " -O3 " CACHE STRING
  " Flags used by the compiler during release builds. ")
set(CMAKE_Fortran_FLAGS_NATIVE " -O4 -march=native -mtune=native -Ofast " CACHE STRING
  " Flags used by the compiler during release builds. ")

IF (TARGET debug)
else()
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug")

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release")

ADD_CUSTOM_TARGET(native
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Native")
endif()
  
# Define runtime output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/EXECUTABLE)
set(SMALL_SFEMaNS_LIB "$ENV{SMALL_SFEMaNS_DIR}/BUILD/libSMALL_SFEMaNS.so")
set(PROJECT_DIR "$ENV{SHALLOW_WATER_DIR}")

file(GLOB_RECURSE global_sources ${PROJECT_DIR}/SOURCES/*.f90)
file(GLOB_RECURSE local_sources ${CMAKE_CURRENT_SOURCE_DIR}/LOCAL_SOURCES/*.f90)

# add executable
set(exe "idp_paraboloid.exe")
add_executable(${exe} ${CMAKE_CURRENT_SOURCE_DIR}/LOCAL_SOURCES/main.f90 ${global_sources} ${local_sources})

# find Intel MKL
set(BLA_VENDOR Intel10_64lp_seq)
find_package(LAPACK REQUIRED)

# include and link library
target_include_directories(${exe} PUBLIC $ENV{MKLROOT}/include $ENV{SMALL_SFEMaNS_DIR}/BUILD)
target_link_libraries(${exe} ${SMALL_SFEMaNS_LIB} LAPACK::LAPACK)

# Regression tests
set(passRegex "ok")

list(APPEND name_of_test    0      1 )
list(APPEND name_of_command low.sh galerkin.sh)

foreach(num IN ZIP_LISTS name_of_test name_of_command)
  set (nm "idp_1D_smooth${num_0}")
  add_test(NAME ${nm}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../REGRESSION_TESTS
  COMMAND /bin/sh/ ${num_1})
  set_property (TEST ${nm} PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
endforeach()
