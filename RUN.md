# Running instructions
Here, we provide a few examples to run and test SFT. Before running the following examples, it is recommended to run the unittest [tests](https://github.com/NOAA-OWP/SoilFreezeThaw/tree/ajk/doc_update/tests).

### Synthetic example (standalone mode)
  - the example uses prescribed soil moisture conditions (static) and ground surface temperature for Laramie, WY. It compares results against the benchmark (golden test results)
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/run_sft.sh">./run_sft.sh</a> STANDALONE (from SoilFreezeThaw directory)    
</pre>

### Real field example (pseudo framework mode)
  - the example runs SFT coupled with Conceptual Funational Equivalent (CFE) and Soil Moisture Profiles (SMP) for about 3 years using Laramie, WY forcing data. The simulated ice_fraction is compared with the existing `golden test` ice_fraction using Schaake runoff scheme. If the test is successful, the user should be able to see `Test passed? Yes`
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/run_sft.sh">./run_sft.sh</a> PFRAMEWORK (from SoilFreezeThaw directory)
</pre>

### Examples (nextgen framework mode)
We provide three examples here to run and test SFT (coupled/uncoupled modes) in the nextgen framework. The following pre-process step needs to be completed before running the below examples. Note run these examples from the ngen directory you cloned [here](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/INSTALL.md#build-nextgen-framework-mode).
  - #### Pre-process step
  ```
  mkdir sft && cd sft
  ln -s ../extern
  ln -s ../data 
  ln -s ./extern/SoilFreezeThaw/SoilFreezeThaw/configs
  ln -s ./extern/SoilFreezeThaw/SoilFreezeThaw/forcings
  ```
  - Example 1: SFT standalone example (uses prescribed ground surface temperature; same as Synthetic test)
  ```
  ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 configs/realization_standalone.json
  ```
  - Example 2: Integrated models (SLoTH+PET+SMP+SFT+CFE) example without NOAH-OWP-Mod (uses prescribed ground surface temperature)
  ```
  ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 configs/realization_multi.json
  ```
  - Example 3: Integrated models (SLoTH+SMP+SFT+CFE/LASAM) example with NOAH-OWP-Mod see [examples](https://github.com/NOAA-OWP/SoilFreezeThaw/tree/ajk/doc_update/examples)