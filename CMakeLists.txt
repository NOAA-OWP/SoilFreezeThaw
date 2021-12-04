cmake_minimum_required(VERSION 3.10)

#set(CMAKE_C_COMPILER $CC)
#set(CMAKE_CXX_COMPILER $CXX)


# set the project name
project(cfeft)

set(CMAKE_BUILD_TYPE Debug)
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)
message(CMAKE_CXX_COMPILER " ${CMAKE_CXX_COMPILER}")
message(CMAKE_C_COMPILER " ${CMAKE_CC_COMPILER}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

# add the executable

## cfe + aorc + ftm
#add_executable(cfeft ./src/main_pass_forcings.cxx ./src/cfe.c ./src/bmi_cfe.c ./forcing_code/src/aorc.c ./forcing_code/src/bmi_aorc.c)

#add_library(cfelib ./externalmodels/bmi_freezethaw/src/bmi_freezethaw.cxx ./externalmodels/bmi_freezethaw/src/freezethaw.cxx ./externalmodels/bmi_freezethaw/include/bmi_freezethaw.hxx ./externalmodels/bmi_freezethaw/include/freezethaw.hxx)

## cfe + aorc + pet + ftm 
add_executable(cfeft ./src/main_cfe_aorc_pet_ftm.cxx ./src/cfe.c ./src/bmi_cfe.c ./forcing_code/src/aorc.c ./forcing_code/src/bmi_aorc.c ./forcing_code/src/pet.c ./forcing_code/src/bmi_pet.c)

add_library(cfelib ./externalmodels/bmi_freezethaw/src/bmi_freezethaw.cxx ./externalmodels/bmi_freezethaw/src/freezethaw.cxx ./externalmodels/bmi_freezethaw/include/bmi_freezethaw.hxx ./externalmodels/bmi_freezethaw/include/freezethaw.hxx)
add_library(freezethaw SHARED ./externalmodels/bmi_freezethaw/src/bmi_freezethaw.cxx ./externalmodels/bmi_freezethaw/src/freezethaw.cxx ./externalmodels/bmi_freezethaw/include/bmi_freezethaw.hxx ./externalmodels/bmi_freezethaw/include/freezethaw.hxx)
add_library(cfebmi SHARED src/bmi_cfe.c src/cfe.c)

target_link_libraries(cfeft LINK_PUBLIC cfelib)
target_link_libraries(cfeft PRIVATE m)


##cmake -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC -DCMAKE_INSTALL_PREFIX=`pwd` -DCMAKE_BUILD_TYPE=Debug ../
##cmake --build .
