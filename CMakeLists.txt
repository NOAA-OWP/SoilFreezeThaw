cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})

option(NGEN "NGEN" OFF)
option(STANDALONE "STANDALONE" OFF)

IF(NGEN)
  message("ngen framework build")
  add_definitions(-DNGEN)
ENDIF(NGEN)

IF(STANDALONE)
  message("Standalone build")
ENDIF(STANDALONE)

# set the project name
project(sftbmi VERSION 1.0.0 DESCRIPTION "OWP SFT BMI Module Shared Library")

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

IF((NOT ${NGEN}) AND (NOT ${STANDALONE}))
  message("Pseudo-framework build")

  add_executable(cfesft ./src/main_cfe_aorc_pet_ftm_smc.cxx ./cfe/src/cfe.c ./cfe/src/giuh.c ./cfe/src/bmi_cfe.c ./cfe/forcing_code/src/aorc.c ./cfe/forcing_code/src/bmi_aorc.c ./pet/src/pet.c ./pet/src/bmi_pet.c)

  add_library(cfesftlib ./src/bmi_soil_freeze_thaw.cxx ./src/soil_freeze_thaw.cxx ./include/bmi_soil_freeze_thaw.hxx ./include/soil_freeze_thaw.hxx ./smc_coupler/src/bmi_soil_moisture_profile.cxx ./smc_coupler/src/soil_moisture_profile.cxx ./smc_coupler/include/bmi_soil_moisture_profile.hxx ./smc_coupler/include/soil_moisture_profile.hxx)

  target_link_libraries(cfesft LINK_PUBLIC cfesftlib)
  target_link_libraries(cfesft PRIVATE m)
ENDIF((NOT ${NGEN}) AND (NOT ${STANDALONE}))

IF(${STANDALONE})
  add_executable(sft ./src/main_sft.cxx ./src/bmi_soil_freeze_thaw.cxx ./src/soil_freeze_thaw.cxx)
ENDIF(${STANDALONE})

# ngen SFT (add shared library)
set(SFT_LIB_NAME_CMAKE sftbmi)
set(SFT_LIB_DESC_CMAKE "OWP SFT BMI Module Shared Library")

# Make sure these are compiled with this directive
add_compile_definitions(BMI_ACTIVE)

if(WIN32)
    add_library(sftbmi src/bmi_soil_freeze_thaw.cxx src/soil_freeze_thaw.cxx)
else()
    add_library(sftbmi SHARED src/bmi_soil_freeze_thaw.cxx src/soil_freeze_thaw.cxx)
endif()

target_include_directories(sftbmi PRIVATE include)

set_target_properties(sftbmi PROPERTIES VERSION ${PROJECT_VERSION})

set_target_properties(sftbmi PROPERTIES PUBLIC_HEADER ./include/bmi_soil_freeze_thaw.hxx)

include(GNUInstallDirs)

install(TARGETS sftbmi
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(sftbmi.pc.in sftbmi.pc @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/sftbmi.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
