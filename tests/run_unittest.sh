#!/bin/bash
set -e

EWTS_PREFIX="${EWTS_PREFIX:-/tmp/ewts_install}"
EWTS_LIB_DIR="${EWTS_LIB_DIR:-${EWTS_PREFIX}/lib}"

if [ ! -d "${EWTS_LIB_DIR}" ] && [ -d "${EWTS_PREFIX}/lib64" ]; then
  EWTS_LIB_DIR="${EWTS_PREFIX}/lib64"
fi

${CXX:-g++} -std=c++17 -Wall -O -g \
  ./main_unittest.cxx \
  ../src/bmi_soil_freeze_thaw.cxx \
  ../src/soil_freeze_thaw.cxx \
  -I../include \
  -I../bmi \
  -I"${EWTS_PREFIX}/include" \
  -L"${EWTS_LIB_DIR}" \
  -Wl,-rpath,"${EWTS_LIB_DIR}" \
  -lewts_cpp \
  -lboost_serialization \
  -lm \
  -o run_sft

./run_sft configs/unittest.txt
test_pass=$?

rm -f run_sft
rm -rf run_sft.dSYM

exit $test_pass
