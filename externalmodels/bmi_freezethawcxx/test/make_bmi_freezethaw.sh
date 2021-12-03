#!/bin/bash
g++ -lm -Wall -O -g ./main_unit_test_bmi.cxx ./../src/bmi_freezethaw.cxx ./../src/freezethaw.cxx -o run_ftm_bmi_test
./run_ftm_bmi_test ./configs/laramie_bmi_config_ftm.txt
