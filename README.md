# Introduction of Soil Freeze-thaw model

## Building the code to run/test examples within the pseudo-framework (not the ngen framework!!)
## Steps to build soil freeze-thaw (SFT) model coupled with cfe and soil moisture profiles
### Setting up (cloning) external repos (cfe and soil moisture profiles repos)
- git clone https://github.com/NOAA-OWP/SoilFreezeThaw && cd SoilFreezeThaw
- git clone https://github.com/NOAA-OWP/cfe 
- git checkout [ajk/cfe_icefraction](https://github.com/NOAA-OWP/cfe/tree/ajk/cfe_icefraction) (note: for the use of SFT in the pseudo-framework we currently use this special branch)
- git clone https://github.com/NOAA-OWP/SoilMoistureProfiles smc_coupler
- mkdir build && build
- cmake -DCMAKE_INSTALL_PREFIX=\`pwd\` -DCMAKE_BUILD_TYPE=Debug ../
- make (to get the executable)
- cd ..
- [run_framework.sh](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/sft_only/run_framework.sh)

# cont.
