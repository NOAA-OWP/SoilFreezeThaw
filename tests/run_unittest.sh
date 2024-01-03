#!/bin/bash
${CXX} -std=c++11 -lm -Wall -O -g ./main_unittest.cxx ../src/bmi_soil_freeze_thaw.cxx ../src/soil_freeze_thaw.cxx -o run_sft
./run_sft configs/unittest.txt
test_pass=$?
rm -f run_sft
rm -rf run_sft.dSYM
exit $test_pass
