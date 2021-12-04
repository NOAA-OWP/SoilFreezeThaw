#!/bin/bash
#CXX --> g++-11 (Homebrew GCC 11.2.0) 11.2.0

##g++ -lm -Wall -O -g ./main_unit_test_bmi.cxx ./../src/bmi_freezethaw.cxx ./../src/freezethaw.cxx -o run_ftm_bmi_test
$CXX -lm -Wall -O0 -g ./main_unit_test_bmi.cxx ./../src/bmi_freezethaw.cxx ./../src/freezethaw.cxx -o run_ftm_bmi_test
./run_ftm_bmi_test ./configs/laramie_bmi_config_ftm.txt
