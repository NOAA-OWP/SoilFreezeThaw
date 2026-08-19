#!/bin/bash
set -e

${CXX:-g++} -std=c++17 -Wall -O -g \
  ./main_unittest.cxx \
  ../src/bmi_soil_freeze_thaw.cxx \
  ../src/soil_freeze_thaw.cxx \
  ../src/Logger.cpp \
  -I../include \
  -I../bmi \
  -lboost_serialization \
  -lm \
  -o run_sft

./run_sft configs/unittest.txt
test_pass=$?

rm -f run_sft
rm -rf run_sft.dSYM

exit $test_pass
