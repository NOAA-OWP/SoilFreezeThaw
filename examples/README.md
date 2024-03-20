# Nextgen SFT Examples
The directory contains two example realization files using CFE and LASAM as soil reservoir models integrated with NOAH-OWP-Mod (NOM), Soil Freeze-thaw (SFT), and Soil Moisture Profiles (SMP). The realizations use the HUC01 catchment (cat-20521) as an example. Follow the instructions below to run these examples. To build the nextgen framework and models see [INSTALL.md](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/INSTALL.md).

**Note:** Your working directory could be anywhere to run these examples, however, it is recommended to work and create the `test_sft` directory out of the `ngen` directory. Assuming you have cloned the SoilFreezeThaw model already, follow the below steps. Also, don't forget to edit `path_to_ngen_dir` lines to adjust them to your local paths.
```
mkdir test_sft
cp -r /path_to_ngen_dir/extern/SoilFreezeThaw/SoilFreezeThaw/examples/ test_sft
cd test_sft
ln -s path_to_ngen_dir ngen
./ngen/cmake_build/ngen ./ngen/data/catchment_data.geojson cat-27 ./ngen/data/nexus_data.geojson nex-26 realizations/realization_<*>.json (where <*> = cfe or lasam)
```