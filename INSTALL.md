# Installation instructions

Detailed instructions on how to build SFT in three modes (standalone, pseudo, and nextgen frameworks). Building SFT requires [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine.

### Build (standalone mode)
 - git clone https://github.com/NOAA-OWP/SoilFreezeThaw
 - cd SoilFreezeThaw && mkdir build && cd build
 - cmake ../ -DSTANDALONE=ON
 - make && cd ..


### Build (pseudo framework mode)
#### Integrated models: [CFE](https://github.com/NOAA-OWP/cfe/), [PET](https://github.com/NOAA-OWP/evapotranspiration), [SMP]( https://github.com/NOAA-OWP/SoilMoistureProfiles), SFT (and more models if needed/desired)
#### Notation: PFRAMEWORK denotes pseudo-framework

 - git clone https://github.com/NOAA-OWP/SoilFreezeThaw && cd SoilFreezeThaw
 - git clone https://github.com/NOAA-OWP/cfe extern/cfe
 - git clone https://github.com/NOAA-OWP/SoilMoistureProfiles extern/SoilMoistureProfiles
 - git clone https://github.com/NOAA-OWP/aorc_bmi extern/aorc_bmi
 - git clone https://github.com/NOAA-OWP/evapotranspiration extern/evapotranspiration
 - mkdir build && cd build
 - cmake ../ -DPFRAMEWORK=ON
 - make && cd ..


### Build (nextgen framework mode)
- See general [instructions](https://github.com/NOAA-OWP/ngen/wiki/NGen-Tutorial#running-cfe) for building models in the ngen framework. 
- ### Specific instructions for building an integrated system
  - git clone https://github.com/noaa-owp/ngen && cd ngen
  - git submodule update --init --recursive
  - #### CFE
    - git submodule update --remote extern/cfe/cfe 
    - cmake -B extern/cfe/cmake_build -S extern/cfe/cfe/ -DNGEN=ON
    - make -C extern/cfe/cmake_build
  - #### fortran bmi
    - cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
    - make -C extern/iso_c_fortran_bmi/cmake_build
  - #### NOAH-OWP-Modular
    - cmake -B extern/noah-owp-modular/cmake_build -S extern/noah-owp-modular
    - make -C extern/noah-owp-modular/cmake_build
  - #### PET
    - cmake -B extern/evapotranspiration/cmake_build -S extern/evapotranspiration/evapotranspiration/
    - make -C extern/evapotranspiration/cmake_build/
  - #### SFT
    - git submodule update --remote extern/SoilFreezeThaw/SoilFreezeThaw  
    - cmake -B extern/SoilFreezeThaw/cmake_build -S extern/SoilFreezeThaw/SoilFreezeThaw/ -DNGEN=ON
    - make -C extern/SoilFreezeThaw/cmake_build
  - #### SMP
    - git submodule update --remote extern/SoilMoistureProfiles/SoilMoistureProfiles
    - cmake -B extern/SoilMoistureProfiles/cmake_build -S extern/SoilMoistureProfiles/SoilMoistureProfiles/ -DNGEN=ON
    - make -C extern/SoilMoistureProfiles/cmake_build
  - cmake -B cmake_build -S . -DBMI_C_LIB_ACTIVE=ON -DBMI_FORTRAN_ACTIVE=ON
  - make -j4 -C cmake_build
  - #### SLoTH is also needed to run SFT in the ngen framework. SLoTH is a BMI that is used to set a bmi variable(s) that is not provided by other BMIs but required by the model. So build [SLoTH](https://github.com/NOAA-OWP/SLoTH) using the following instructions
    - cd extern/sloth/ && git checkout latest 
    - git submodule update --init --recursive
    - cd ../..
    - cmake -B extern/sloth/cmake_build -S extern/sloth/
    - make -C extern/sloth/cmake_build
