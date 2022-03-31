cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})

option(NGEN "NGEN" OFF)
IF(NGEN)
message("ngen framework build")
ENDIF(NGEN)


# set the project name
project(cfesft)

set(CMAKE_BUILD_TYPE Debug)
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)
message(CMAKE_CXX_COMPILER " ${CMAKE_CXX_COMPILER}")
message(CMAKE_C_COMPILER " ${CMAKE_C_COMPILER}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

# add the executable

## cfe + aorc + pet + ftm

IF(NOT ${NGEN})
message("Pseudo-framework build")
add_executable(cfesft ./src/main_cfe_aorc_pet_ftm_smc.cxx ./cfe/src/cfe.c ./cfe/src/bmi_cfe.c ./cfe/forcing_code/src/aorc.c ./cfe/forcing_code/src/bmi_aorc.c ./cfe/forcing_code/src/pet.c ./cfe/forcing_code/src/bmi_pet.c)
add_library(cfesftlib ./src/bmi_freezethaw.cxx ./src/freezethaw.cxx ./include/bmi_freezethaw.hxx ./include/freezethaw.hxx ./smc_coupler/src/bmi_coupler.cxx ./smc_coupler/src/smc_profile.cxx ./smc_coupler/include/bmi_coupler.hxx ./smc_coupler/include/smc_profile.hxx)
target_link_libraries(cfesft LINK_PUBLIC cfesftlib)
target_link_libraries(cfesft PRIVATE m)
ENDIF(NOT ${NGEN})


# ngen (add shared library)
add_library(sftbmi SHARED ./src/bmi_freezethaw.cxx ./src/freezethaw.cxx ./include/bmi_freezethaw.hxx ./include/freezethaw.hxx)


target_include_directories(sftbmi PRIVATE /include)

##cmake -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC -DCMAKE_INSTALL_PREFIX=`pwd` -DCMAKE_BUILD_TYPE=Debug ../
##cmake --build .



