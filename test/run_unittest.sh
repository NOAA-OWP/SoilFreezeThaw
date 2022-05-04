#!/bin/bash
${CXX} -lm -Wall -O -g ./main_unit_test_bmi.cxx ../src/bmi_soil_freeze_thaw.cxx ../src/soil_freeze_thaw.cxx -o run_sft
./run_sft configs/unittest.txt
