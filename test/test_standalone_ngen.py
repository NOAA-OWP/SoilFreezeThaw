import numpy as np
import pandas as pd

golden_file="extern/SoilFreezeThaw/SoilFreezeThaw/test/file_golden.csv"

ngen_file="cat-27.csv"

total_num_hours = 24566


dat_ngen_icef = pd.read_csv(ngen_file,sep=',',skiprows=0, nrows=total_num_hours)

dat_standalone_icef = pd.read_csv(golden_file,sep=',',skiprows=0, nrows=total_num_hours)

err_frozen_frac_mm = 0
for i in range(total_num_hours):
    errorerr_frozen_frac_mm = err_frozen_frac_mm + dat_ngen_icef["ice_fraction_schaake"] - dat_standalone_icef["ice_fraction"]

#RMSE
err_frozen_frac_mm  = pow(err_frozen_frac_mm/total_num_hours,0.5);

test_status = False
if (err_frozen_frac_mm < 1.e-2):
    test_status = True
else:
    test_status = False

passed = "YES" if test_status else "No"


print("*********************************************************\n")
print("    Test passed = ", passed)
print("    Ice fraction error = ",err_frozen_frac_mm)
print("*********************************************************")
