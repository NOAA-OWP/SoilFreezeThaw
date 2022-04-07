# Unit test for Soil Freeze Thaw model
Usage: run `./make_bmi_freezethaw.sh`

Multiple checks are performed:
1. Check number of input/output variables
2. Test `GetVar*` methods for the input variables and compare against initial data; if failed, will throw an error
3. Test `GetGrid*` methods and compare data against known values, such as endtime etc.
4. Loop over the input variables, use `Set*` and `Get*` methods to verify `Get*` return the same data set by `Set*`
5. Step (4) for output variables
6. Using `Update` method, take 3333 timestep and compare the `soil temperature` and `frozen fraction` against benchmark test