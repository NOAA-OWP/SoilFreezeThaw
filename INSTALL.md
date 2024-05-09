# Build and Run Instructions
Detailed instructions on how to build and run SoilFreezeThaw (SFT) in three modes (standalone, pseudo, and nextgen frameworks) are provided below. Standalone mode refers to running the model by itself, not linked to any other model.  Pseudo or pseudo framework mode uses code within the SoilFreezeThaw repo to link external models with SFT.  The nextgen framework mode creates a build that must be run within the [nextgen framework](https://github.com/NOAA-OWP/ngen), a framework used to ingest data, pass data, link models, and provide model output.

Building SFT requires [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine.

## Clone the repository
```
git clone https://github.com/NOAA-OWP/SoilFreezeThaw
cd SoilFreezeThaw 
```
***Note:** Before running the following examples, it is recommended to run the unittests [tests](https://github.com/NOAA-OWP/SoilFreezeThaw/tree/master/tests).

## Standalone mode example
The standalone example uses prescribed soil moisture conditions that do not change with time (static) and ground surface temperature for Laramie, WY.  It is intended for testing purposes only. The simulated results are compared against the benchmark results.

### Build
```
 mkdir build && cd build
 cmake ../ -DSTANDALONE=ON
 make && cd ..
```
### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/run_sft.sh">./run_sft.sh</a> STANDALONE (from SoilFreezeThaw directory)    
</pre>

## Pseudo framework mode example
The pseudo framework example runs SFT coupled with Conceptual Funational Equivalent [CFE](https://github.com/NOAA-OWP/cfe/), Soil Moisture Profiles [SMP]( https://github.com/NOAA-OWP/SoilMoistureProfiles), and potential evapotranspiration model [PET](https://github.com/NOAA-OWP/evapotranspiration) for about 3 years using Laramie, WY [forcing data](https://github.com/NOAA-OWP/SoilFreezeThaw/tree/master/forcings). 

The simulated ice_fraction is compared with the existing `golden test` ice_fraction using Schaake runoff scheme. If the test is successful, the user should be able to see `Test passed? Yes`.
**Notation:*** PFRAMEWORK denotes pseudo-framework
### Build
All these steps should happen inside the SoilFreezeThaw directory.
```
git clone https://github.com/NOAA-OWP/cfe extern/cfe
git clone https://github.com/NOAA-OWP/SoilMoistureProfiles extern/SoilMoistureProfiles
git clone https://github.com/NOAA-OWP/aorc_bmi extern/aorc_bmi
git clone https://github.com/NOAA-OWP/evapotranspiration extern/evapotranspiration
mkdir build && cd build
cmake ../ -DPFRAMEWORK=ON
make && cd ..
```
### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/run_sft.sh">./run_sft.sh</a> PFRAMEWORK (from SoilFreezeThaw directory)
</pre>

## Nextgen framework mode examples
We provide three examples here to run and test SFT (coupled/uncoupled modes) in the nextgen framework.
- Example 1: SFT standalone example (uses prescribed ground surface temperature; same as Synthetic test). This method should only be used for testing purposes.
- Example 2: Integrated models (SLoTH+PET+SMP+SFT+CFE) example without NOAH-OWP-Mod (uses prescribed ground surface temperature). 
- Example 3: Integrated models (SLoTH+SMP+SFT+CFE/LASAM) example with NOAH-OWP-Mod where NOAH-OWP-Mod provides the ground surface temperature to SFT.

### Build
See general [instructions](https://github.com/NOAA-OWP/ngen/wiki/NGen-Tutorial#running-cfe) for building models in the ngen framework. Note that ngen cmake flag can be set (i.e., `NGEN_WITH_EXTERN_ALL=ON`) to build SLoTH, CFE, PET, Noah-owp-modular (more details are provided below).

Specific instructions for building an integrated system:

#### 1. Clone ngen & extract submodules
```
git clone https://github.com/noaa-owp/ngen && cd ngen
git submodule update --init --recursive
```
#### 2. Fortran BMI
To use the noah-owp-modular fortran library, we need to build the iso_c_fortran_bmi middleware.
```
cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
make -C extern/iso_c_fortran_bmi/cmake_build
```
#### 3. CFE
```
git submodule update --remote extern/cfe/cfe 
cmake -B extern/cfe/cfe/cmake_build -S extern/cfe/cfe/ -DNGEN=ON
make -C extern/cfe/cfe/cmake_build
```
#### 4. LASAM
```
git clone https://github.com/NOAA-OWP/LGAR-C extern/LASAM
cmake -B extern/LASAM/LASAM/cmake_build -S extern/LASAM/ -DNGEN=ON
make -C extern/LASAM/LASAM/cmake_build
```
#### 5. NOAH-OWP-Modular
```
cmake -B extern/noah-owp-modular/noah-owp-modular/cmake_build -S extern/noah-owp-modular -DNGEN_IS_MAIN_PROJECT=ON
make -C extern/noah-owp-modular/noah-owp-modular/cmake_build
```
#### 6. PET
```
cmake -B extern/evapotranspiration/evapotranspiration/cmake_build -S extern/evapotranspiration/evapotranspiration/
make -C extern/evapotranspiration/evapotranspiration/cmake_build/
```
#### 7. SFT
```
git submodule update --remote extern/SoilFreezeThaw/SoilFreezeThaw  
cmake -B extern/SoilFreezeThaw/SoilFreezeThaw/cmake_build -S extern/SoilFreezeThaw/SoilFreezeThaw/ -DNGEN=ON
make -C extern/SoilFreezeThaw/SoilFreezeThaw/cmake_build
```
#### 8. SMP
```
git submodule update --remote extern/SoilMoistureProfiles/SoilMoistureProfiles
cmake -B extern/SoilMoistureProfiles/SoilMoistureProfiles/cmake_build -S extern/SoilMoistureProfiles/SoilMoistureProfiles/ -DNGEN=ON
make -C extern/SoilMoistureProfiles/SoilMoistureProfiles/cmake_build
```
#### 9. SLoTH
SLoTH is also needed to run SFT in the ngen framework. SLoTH is a BMI that is used to set a bmi variable(s) that is not provided by other BMIs but required by the model. So build [SLoTH](https://github.com/NOAA-OWP/SLoTH) using the following instructions.
```
cd extern/sloth/ && git checkout latest 
git submodule update --init --recursive
cd ../..
cmake -B extern/sloth/cmake_build -S extern/sloth/
make -C extern/sloth/cmake_build
```
#### 10. Build ngen
```
cmake -B cmake_build -S . -DNGEN_WITH_BMI_C=ON -DNGEN_WITH_BMI_FORTRAN=ON -DNGEN_WITH_SQLITE=ON
make -j4 -C cmake_build
```

As mentioned earlier, CFE, SLoTH, PET, and Noah-owp-modular can be built as ngen external models at the ngen build time using
```
cmake -B cmake_build -S . -DNGEN_WITH_BMI_FORTRAN=ON -DNGEN_WITH_SQLITE=ON -DNGEN_WITH_EXTERN_ALL=ON
make -j4 -C cmake_build
```

### Run
The following pre-process step needs to be completed before running the below examples.
  #### Pre-process step
  ```
  mkdir sft && cd sft
  ln -s ../extern
  ln -s ../data 
  ln -s ./extern/SoilFreezeThaw/SoilFreezeThaw/configs
  ln -s ./extern/SoilFreezeThaw/SoilFreezeThaw/forcings
  ln -s ./extern/SoilFreezeThaw/SoilFreezeThaw/realizations
  ```
  
  #### Example 1
  ```
  ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realizations/realization_standalone.json
  ```
  Run: `python ./extern/SoilFreezeThaw/SoilFreezeThaw/tests/test_standalone_ngen.py` (from the sft directory) to compare Example 1 results against the benchmark.
  #### Example 2
  ```
  ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realizations/realization_coupled.json
  ```
  #### Example 3
  Detailed instructions for running this example along with realization and model config files are provided in the [examples](https://github.com/NOAA-OWP/SoilFreezeThaw/tree/master/examples) directory.
## Post-process step
  - For standalone simulations: run `python extern/SoilFreezeThaw/SoilFreezeThaw/tests/test_standalone_ngen.py` ([test_standalone_ngen.py](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/tests/test_standalone_ngen.py) script compares results with a gold test [output](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/tests/file_golden.csv))
  - For integrated simulation: Output data is stored in cat-27.csv, use your favorite tool to visualize data
