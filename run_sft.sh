#!/bin/bash
#g++ -lm -Wall -O -g ./src/main_sft.cxx ./src/bmi_freezethaw.cxx ./src/freezethaw.cxx -o run_sft
./build/sft configs/laramie_bmi_config_ftm_standalone.txt
