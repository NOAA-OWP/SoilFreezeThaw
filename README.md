## Standalone soil freeze-thaw (SFT) model example
### Example description: 
Runs SFT for about 3 years using Laramie, WY forcing data. The simulated ice_fraction is compared with existing `golden test` ice_fraction using Schaake scheme. If test is successfull, the user should be able to see `Test passed = Yes` 
### Build
- git clone https://github.com/NOAA-OWP/SoilFreezeThaw && cd SoilFreezeThaw
- mkdir build && cd build
- cmake -DCMAKE_INSTALL_PREFIX=\`pwd\` -DCMAKE_BUILD_TYPE=Debug -DSTANDALONE=ON ../
- make && cd ..
### Run
- [run_sft.sh](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/run_sft.sh)
- It compares results with gold test result


## Pseudo-framework integrated models example 
### Integrated models: [CFE](https://github.com/NOAA-OWP/cfe/), [PET](https://github.com/NOAA-OWP/evapotranspiration), [SMP]( https://github.com/NOAA-OWP/SoilMoistureProfiles), SFT (and more models if needed/desired)
### Build 
- git clone https://github.com/NOAA-OWP/SoilFreezeThaw && cd SoilFreezeThaw
- git clone https://github.com/NOAA-OWP/cfe 
- git checkout [cfe_soilfreezethaw](https://github.com/NOAA-OWP/cfe/tree/cfe_soilfreezethaw) (note: coupling SFT with CFE using the pseudo-framework requires this special cfe branch)
- git clone https://github.com/NOAA-OWP/SoilMoistureProfiles smc_coupler
- git clone https://github.com/NOAA-OWP/evapotranspiration pet
- mkdir build && cd build
- cmake -DCMAKE_INSTALL_PREFIX=\`pwd\` -DCMAKE_BUILD_TYPE=Debug ../
- make && cd ..
### Run
- [run_framework.sh](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/run_framework.sh)

## ngen-framework standalone/integrated models example
- See general [instructions](https://github.com/NOAA-OWP/ngen/wiki/NGen-Tutorial#running-cfe) for building models in the ngen framework. 
- ### Specific instructions for building an integrated system
### Build
  - git clone https://github.com/noaa-owp/ngen && cd ngen
  - git submodule update --init
  - #### CFE
    - **Note:** make sure to pull latest cfe and checkout cfe_soilfreezethaw branch: 
    - git submodule update --remote extern/cfe/cfe  
    - cd extern/cfe/cfe && git checkout cfe_soilfreezethaw && cd ../../..
    - cmake -B extern/cfe/cmake_build -S extern/cfe
    - make -C extern/cfe/cmake_build
  - #### fortran bmi
    - cmake -B extern/iso_c_fortran_bmi/cmake_build -S extern/iso_c_fortran_bmi
    - make -C extern/iso_c_fortran_bmi/cmake_build
  - #### PET
    - cmake -B extern/evapotranspiration/cmake_build -S extern/evapotranspiration/evapotranspiration/
    - make -C extern/evapotranspiration/cmake_build/
  - #### SFT
    - git submodule update --remote extern/SoilFreezeThaw/SoilFreezeThaw  
    - cmake -B extern/SoilFreezeThaw/cmake_build -S extern/SoilFreezeThaw/SoilFreezeThaw/ -DNGEN:BOOL=ON
    - make -C extern/SoilFreezeThaw/cmake_build
  - #### SMP
    - git submodule update --remote extern/SoilMoistureProfiles/SoilMoistureProfiles
    - cmake -B extern/SoilMoistureProfiles/cmake_build -S extern/SoilMoistureProfiles/SoilMoistureProfiles/ -DNGEN:BOOL=ON
    - make -C extern/SoilMoistureProfiles/cmake_build 
  - Note the following step will be removed once [PR#405](https://github.com/NOAA-OWP/ngen/pull/405) is merged
    - gh pr list (you should see #405  Bmi array support  hellkite500:bmi-array-support)
    - gh pr checkout 405
  - cmake -B cmake_build -S . -DNGEN_ACTIVATE_PYTHON:BOOL=ON -DBMI_C_LIB_ACTIVE:BOOL=ON -DBMI_FORTRAN_ACTIVE:BOOL=ON
  - make -j4 -C cmake_build

 ### Run 
 #### Pre-process step
 ```
   mkdir sft && cd sft
   ln -s ../extern
   ln -s ../data 
  ```
 #### standalone SFTM in the ngen framework
 ```
   cp extern/SoilFreezeThaw/SoilFreezeThaw/configs/realization_config_standalone_<macos/linux>.json . (pick the file based on your machine)
   ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realization_config_standalone.json
 ```
#### Run integrated models in the ngen framework
 ```
   cp extern/SoilFreezeThaw/SoilFreezeThaw/configs/realization_config_multi_<macos/linux>.json . (pick the file based on your machine)
   ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realization_config_multi_<macos/linux>.json
  ```
#### Post-process step
  - For standalone simulaition: run `python extern/SoilFreezeThaw/SoilFreezeThaw/test/test_standalone_ngen.py` ([test_standalone_ngen.py](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/test/test_standalone_ngen.py) script compares results with a gold test [output](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/test/file_golden.csv))
  - For integrated simulation: Output data is stored in cat-27.csv, use your favorite tool to visualize data (or to compare with pseudo framework)


## Parameters in the config file
* `forcing_file` provides ground temperatures used as a surface boundary condition. Won't be needed when coupled to models (e.g., noah-om) that can provide ground temperature data
* `smcmax` (porosity; maximum moisture content)
* `b` is the pore size distribution, beta exponent on Clapp-Hornberger
* `satpsi` saturated capillary head (saturated moisture potential)
* `quartz` soil quartz content, used in soil thermal conductivity function of Peters-Lidard
* `soil_z` is the vertical resolution of the soil column (computational domain)
* `soil_temperature` is the initial soil temperature for the discretized column
* `soil_moisture_content` is the initial soil total (liquid + ice) moisture content for the discretized column
* `soil_liquid_content` is the initial soil liquid moisture content for the discretized column
* `ice_fraction_scheme` is associated with the runoff scheme used in the soil reservoir models (e.g. CFE), options: Schaake and Xinanjiang
* `bottom_boundary_temp` sets the temperature conditions at the bottom boundary of the domain/column, if not specified zero-geothermal boundary condition is used
* `sft_standalone` options: true or false. True for soil freeze-thaw model standalone run.
* `soil_moisture_bmi` options: true or false. If true soil_moisture profile is set by the SoilMoisutreProfile module throught its BMI.

## Introduction of Soil Freeze-thaw model

The diffusion equation is used to simulate the transport of energy in the soil.

<img width="615" alt="eq1" src="https://user-images.githubusercontent.com/15165757/157314534-0ef6e5ea-4dad-4be5-aaac-888ca139cbc7.png">

where T is the soil temperature, C is the volumetric heat capacity [J/(m3K)], ?? is the volumetric soil moisture content [m3/m3], K is the bulk thermal conductivity [W/(mK)], Lf is the latent heat of fusion (J/Kg), ice is the ice density, ice is the soil ice content [m3/m3], t is the time, and z is the depth. The last term on the right side of Eq. (1) represents the latent heat of released/consumed during the phase change. The soil thermal conductivity follows the parameterization of Peters-Lidars (Ref) given by:

<img width="598" alt="eq2" src="https://user-images.githubusercontent.com/15165757/157315994-bca4d024-0288-4809-aeae-92e9012b4318.png">

where Ke , Ksat, and Kdry are the Kersten number, dry, and saturated soil thermal conductivities, respectively, and are defined in (Ref). The effective volumetric heat capacity, C, is calculated based on the respective fraction of each component (water, ice, air, and rock):

<img width="617" alt="eq3" src="https://user-images.githubusercontent.com/15165757/157316057-3aeffa5e-1b0e-4a29-b34a-12bda5d1e58e.png">

here, ??, ??_w, and ??_s are the total water content, the liquid water content and the maximum soil moisture content (porosity), respectively. The parameters  Cw, Cice, Cair, and Csoil are the volumetric heat capacities of water, ice, air, and soil, respectively. Table?
The diffusion equation (Eq. 1) is discretized using the Crank-Nicolson scheme (Ref), which is an implicit finite difference scheme and thus unconditionally stable.

Freezing-point depression equation 
Frozen soils not only impact the soil thermal and hydrological properties, but its presence strongly affects the partition of precipitation into surface runoff and infiltration leading to a strongly coupled surface and subsurface system. Due to the formation of impermeable soil over the winter period, infiltration of water significantly reduces during the spring freshet generating high peaks in the surface discharge.

In unfrozen conditions, the liquid water content is equal to the total water content in the soil. However, when soil freezes the liquid water content also depends on the soil temperature. The soil matric potential,  [m], under the assumption that soil water potential remains in equilibrium with vapor pressure over pure ice when ice is present in the soil, and neglecting soil water osmotic potential, is given by:

<img width="607" alt="eq4" src="https://user-images.githubusercontent.com/15165757/157316409-3c51afda-e64b-4f6b-bf9e-efa50466cae4.png">	

where T and To  are the soil temperature and freezing temperatures, respectively (Ref), g is the acceleration of gravity [m s-2] and Lf is the latent heat of fusion [J Kg-1]. The soil matric potential, psi [m], as a function of soil moisture content is given by (ref: Clapp-Hornberger):

<img width="594" alt="eq5" src="https://user-images.githubusercontent.com/15165757/157317205-a2b1e455-8714-4113-8c56-f7a2ef55533c.png">

where ??_s,??_s,and ??_l are saturated soil matric potential [m], saturated soil moisture [-] (porosity), and soil liquid moisture [-], respectively. The exponent b is the Clapp-Hornberger parameter. According to the ???freezing equals drying??? approximation (Ref), we equate Eq. (4) and Eq. (5) to obtain an expression for liquid water content:

<img width="626" alt="eq6" src="https://user-images.githubusercontent.com/15165757/157317016-3b8c0ad1-1f4b-4370-9773-c97aba0d0ee5.png">

Equation (6) is called the ???freezing-point depression equation??? (Refs) and gives the maximum amount of liquid water (unfrozen soil moisture content) that can exist below the subfreezing temperature. The frozen soil moisture content (i, ice fraction) is given by:
??_ice = ??_t- ??_l , where t is the total water content given by soil water retention curve: 

<img width="385" alt="eq8" src="https://user-images.githubusercontent.com/15165757/158890563-0ff39857-c8d4-41fe-ab06-2f5548ec25a1.png">
