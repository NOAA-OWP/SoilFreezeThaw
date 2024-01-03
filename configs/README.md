## Configuration File
Example configuration files are provided in this directory. To build and run the given examples see the instructions [here](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/INSTALL.md).

A detailed description of the parameters for model configuration (i.e., initialize/setup) is provided below. The asterisk (*) denotes calibratable parameters, (i.e., `smcmax, b, and satpsi` can be calibrated)
| Variable ______________ | Datatype ________ | Limits ______ | Units ______ | Role _____ |  Description __________________________________________________ |
| ------ | -------- | ------ | ----- | ---- | ----------- |
| forcing_file | string | - | - | filename | provides ground temperature (not needed when coupled to models providing ground temperature data|
| *soil_params.smcmax | double | - | - | state variable | maximum soil moisture content (porosity) |
| *soil_params.b | double | - | m | state variable | pore size distribution, beta exponent in Clapp-Hornberger characteristic function |
| *soil_params.satpsi | double | - | m | state variable | saturated capillary head (saturated moisture potential) |
| soil_params.quartz | double | - | m | state variable | soil quartz content, used in soil thermal conductivity function of Peters-Lidard |
| ice_fraction_scheme | int | - | - | coupling variable | runoff scheme used in the soil reservoir models (e.g. CFE), options: Schaake and Xinanjiang|
| soil_z | double (1D array) | - | m | spatial resolution | vertical resolution of the soil column (computational domain of the SFT model) |
| soil_temperature | double (1D array) | - | K | spatial resolution | initial soil temperature for the discretized column |
| soil_moisture_content | double (1D array) | - | - | spatial resolution | initial soil total (liquid + ice) moisture content for the discretized column |
| soil_liquid_content | double (1D array) | - | - | spatial resolution | initial soil liquid moisture content for the discretized column|
| bottom_boundary_temp | double | - | K | boundary condition | temperature at the bottom boundary (BC) of the domain, if not specified, the default BC is zero-geothermal flux|
| top_boundary_temp | double | - | K | boundary condition | temperature at the top/surface boundary of the domain, if not specified, then other options include: 1) read from a file, or 2) provided through coupling |
| sft_standalone | boolean | true, false | - | coupling variable | true for standalone model run; default is false |
| soil_moisture_bmi | boolean | true, false | - | coupling variable | If true soil_moisture_profile is set by the SoilMoistureProfile module through the BMI; if false then config file must provide soil_moisture_content and soil_liquid_content |
| dt | double | - | s, sec, h, hr, d, day | timestep size | Size of a simulation timestep. If no unit is specified defaults to hour. |
| end_time | double | - | s, sec, h, hr, d, day  | simulation duration | Simulation duration. If no unit is specified defaults to hour. |
