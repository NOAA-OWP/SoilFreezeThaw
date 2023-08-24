# Soil Freeze-thaw Model
The soil freeze-thaw model simulates the transport of heat in soil using a one-dimensional vertical column. The model uses a standard diffusion equation discretized using a fully-implicit scheme at the interior and a semi-implicit scheme at the top and bottom boundaries, similar to NOAH-MP. More details are provided below.

## Build and Run Instructions
Detailed instructions on how to build and run SoilFreezeThaw (SFT) model can be found in [INSTALL](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/INSTALL.md) guide.
  - Test examples highlights
    - Unittest (see [tests](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/tests/README.md))
    - Synthetic example (standalone mode): simulations with prescribed soil moisture profiles (static input) (see [build/run](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/INSTALL.md#standalone-mode-example))
    - Real field example (pseudo framework mode): simulations with real forcing data (see [build/run](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/INSTALL.md#pseudo-framework-mode-example))
    - Examples (nextgen framework mode):
       - Synthetic example: Identical to the above `Synthetic example` but runs in the nextgen framework.
       - Real field example: Identical to the above `Real field example` but runs in the nextgen framework.
       - Real field example: Two nextgen realization examples coupling 1) SFT with [CFE](https://github.com/NOAA-OWP/cfe/) and 2) SFT with [LASAM](https://github.com/NOAA-OWP/LGAR-C) running on a catchment in HUC01 region are also provided in [examples](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/examples/README.md) directory. Build and run instructions are given at [build](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/INSTALL.md#nextgen-framework-mode-example) and [run](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/examples/README.md), respectively.
    
## Model Configuration File
A detailed description of the parameters for model configuration is provided [here](https://github.com/NOAA-OWP/SoilFreezeThaw/tree/ajk/doc_update/configs/README.md).
  
## Getting help
For questions, please contact Ahmad Jan (ahmad.jan(at)noaa.gov), the main developer/maintainer of the repository.

## Known issues or raise an issue
We are constantly looking to improve the model and/or fix bugs as they arise. Please see the Git Issues for known issues or if you want to suggest adding a capability or to report a bug, please open an issue.

## Getting involved
See general instructions to contribute to the model development ([instructions](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/CONTRIBUTING.md)). Simply fork the repository and submit a pull request.
