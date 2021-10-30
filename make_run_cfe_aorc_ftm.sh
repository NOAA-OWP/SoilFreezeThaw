#!/bin/bash
##gcc -lm ./src/main_pass_forcings.c ./src/cfe.c ./src/bmi_cfe.c ./forcing_code/src/aorc.c ./forcing_code/src/bmi_aorc.c ./externalmodels/bmi_freezethaw/src/bmi_freezethaw.c ./externalmodels/bmi_freezethaw/src/freezethaw.c -o run_cfe_bmi_pass_forcingsxx
_build/cfeft ./configs/cat_87_bmi_config_cfe_pass.txt ./configs/cat_87_bmi_config_aorc.txt ./configs/config_ftm.txt
