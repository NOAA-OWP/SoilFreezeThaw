cmake_minimum_required(VERSION 3.10)

#set(CMAKE_C_COMPILER $CC)
#set(CMAKE_CXX_COMPILER $CXX)
set(CMAKE_BUILD_TYPE Debug)

# set the project name
project(cfeft)


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_CXX_FLAGS "-std=c++0x")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
##set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libcxx -lc++abi")
# add the executable


## cfe + aorc + ftm
#add_executable(cfeft ./src/main_pass_forcings.cxx ./src/cfe.c ./src/bmi_cfe.c ./forcing_code/src/aorc.c ./forcing_code/src/bmi_aorc.c)

#add_library(cfelib ./externalmodels/bmi_freezethawcxx/src/bmi_freezethaw.cxx ./externalmodels/bmi_freezethawcxx/src/freezethaw.cxx ./externalmodels/bmi_freezethawcxx/include/bmi_freezethaw.hxx ./externalmodels/bmi_freezethawcxx/include/freezethaw.hxx)

## cfe + aorc + pet + ftm 
add_executable(cfeft ./src/main_cfe_aorc_pet_ftm.cxx ./src/cfe.c ./src/bmi_cfe.c ./forcing_code/src/aorc.c ./forcing_code/src/bmi_aorc.c ./forcing_code/src/pet.c ./forcing_code/src/bmi_pet.c)

add_library(cfelib ./externalmodels/bmi_freezethawcxx/src/bmi_freezethaw.cxx ./externalmodels/bmi_freezethawcxx/src/freezethaw.cxx ./externalmodels/bmi_freezethawcxx/include/bmi_freezethaw.hxx ./externalmodels/bmi_freezethawcxx/include/freezethaw.hxx)

target_link_libraries(cfeft LINK_PUBLIC cfelib)
target_link_libraries(cfeft PRIVATE m)


##cmake -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC -DCMAKE_INSTALL_PREFIX=`pwd` -DCMAKE_BUILD_TYPE=Debug ../
##cmake --build .
