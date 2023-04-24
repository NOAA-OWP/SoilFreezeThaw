cmake_minimum_required(VERSION 3.10)

string(ASCII 27 Esc)
set(ColourReset "${Esc}[m")
set(Red         "${Esc}[32m")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})

option(NGEN "NGEN" OFF)
option(STANDALONE "STANDALONE" OFF)
option(STANDALONE "PFRAMEWORK" OFF)

if(NGEN)
  message("ngen framework build")
  add_definitions(-DNGEN)
endif(NGEN)

if(NOT STANDALONE AND NOT PFRAMEWORK AND NOT NGEN)
  message("${Red}Options: STANDALONE, PFRAMEWORK, NGEN" ${ColourReset})
  message(FATAL_ERROR "Invalid option is provided, CMake will exit." )
endif()


# set the project name
project(sftbmi VERSION 1.0.0 DESCRIPTION "OWP SFT BMI Module Shared Library")

set(CMAKE_BUILD_TYPE Debug)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
  message("Debug build.")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

message("CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message("CMAKE_C_COMPILER   = ${CMAKE_C_COMPILER}")
message("CMAKE_BUILD_TYPE   = ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

if(STANDALONE)
message("${Red} Soil freeze-thaw model standalone build! ${ColourReset}")
set(exe_name "sft_standalone")
elseif(PFRAMEWORK)
message("${Red} Soil freeze-thaw model pseudo-framework build! ${ColourReset}")
set(exe_name "sft_pframework")
elseif(NGEN)
message("${Red} Soil freeze-thaw model ngen-framework build! ${ColourReset}")
endif()

# add the executable

## cfe + aorc + pet + ftm

if(PFRAMEWORK)
  add_executable(${exe_name} ./src/main_pseudo_framework.cxx
                 ./extern/cfe/src/cfe.c ./extern/cfe/src/bmi_cfe.c
		 ./extern/cfe/src/giuh.c ./extern/cfe/src/conceptual_reservoir.c
		 ./extern/aorc_bmi/src/aorc.c ./extern/aorc_bmi/src/bmi_aorc.c
		 ./extern/evapotranspiration/src/pet.c ./extern/evapotranspiration/src/bmi_pet.c)

  add_library(sftlib ./src/bmi_soil_freeze_thaw.cxx ./src/soil_freeze_thaw.cxx
              ./include/bmi_soil_freeze_thaw.hxx ./include/soil_freeze_thaw.hxx
	      ./extern/SoilMoistureProfiles/src/bmi_soil_moisture_profile.cxx
	      ./extern/SoilMoistureProfiles/src/soil_moisture_profile.cxx
	      ./extern/SoilMoistureProfiles/include/bmi_soil_moisture_profile.hxx
	      ./extern/SoilMoistureProfiles/include/soil_moisture_profile.hxx)

  target_link_libraries(${exe_name} LINK_PUBLIC sftlib)
  target_link_libraries(${exe_name} PRIVATE m)
  target_include_directories(${exe_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/extern/cfe/include)
  target_include_directories(${exe_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/extern/)
elseif(STANDALONE)
  add_executable(${exe_name} ./src/main_standalone.cxx ./src/bmi_soil_freeze_thaw.cxx ./src/soil_freeze_thaw.cxx)
endif()

##for NGEN BUILD
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


unset(STANDALONE CACHE)
unset(PFRAMEWORK CACHE)
unset(NGEN CACHE)
