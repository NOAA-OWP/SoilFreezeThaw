The directory contains two example realization files using CFE and LASAM as soil revervoir models integrated with NOAH-OWP-Mod (NOM), 
Soil Freeze-thaw (SFT), and Soil Moisture Profiles (SMP). The realizations use HUC01 catchment (cat-20521) as an example. Follow the instructions
below to run these examples.

Build ngen framework and models (see [README.md](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/README.md) for instructions)
```
cp -r /path_to_ngen_dir/extern/SoilFreezeThaw/SoilFreezeThaw/examples test
cd test
ln -s /path_to_ngen_dir/extern
ln -s /path_to_ngen_dir/data
ln -s /path_to_ngen_dir/cmake_build
./cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realizations/realization_X.json (where X = cfe or lasam)
```
