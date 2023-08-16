## Soil Freeze-thaw Model
The soil freeze-thaw model simulates the transport of heat in soil using one-dimensional vertical column. The model uses standard diffusion equation discretized using fully-implicit scheme at the interior and semi-implicit scheme at the top and bottom boundaries, similar to NOAH-MP. More details are provided below.

### Instructions for: 
  - Installation (see [instructions](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/INSTALL.md))
  - Test examples
    - Unittest (see [tests](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/tests/README.md))
    - Synthetic example: simulations with synthetic forcing data (see [run](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/RUN.md#synthetic-example-standalone-mode))
    - Real field example: simulations with real forcing data (see [run](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/RUN.md#real-field-example-pseudo-framework-mode))

### Model Configuration File
  - Detailed description of the parameters for model configuration is provided ([here](https://github.com/NOAA-OWP/SoilFreezeThaw/tree/ajk/doc_update/configs/README.md))

### Nextgen Realization Files
  - Realization files for running LASAM (coupled/uncoupled modes) in the nextgen framework are provided here ([here](https://github.com/NOAA-OWP/SoilFreezeThaw/tree/ajk/doc_update/realizations/README.md))
  
### Getting help
For questions, please contact Ahmad Jan (ahmad.jan(at)noaa.gov), the main developer/maintainer of the repository.

### Known issues or raise an issue
We are constantly looking to improve the model and/or fix bugs as they arise. Please see the Git Issues for known issues or if you want to suggest adding a capability or to report a bug, please open an issue.

### Getting involved
See general instructions to contribute to the model development ([instructions](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/CONTRIBUTING.md)). Simply fork the repository and submit a pull request.

## Introduction of Soil Freeze-thaw model

The diffusion equation is used to simulate the transport of energy in the soil.

<img width="615" alt="eq1" src="https://user-images.githubusercontent.com/15165757/157314534-0ef6e5ea-4dad-4be5-aaac-888ca139cbc7.png">

where T is the soil temperature, C is the volumetric heat capacity [J/(m3K)], θ is the volumetric soil moisture content [m3/m3], K is the bulk thermal conductivity [W/(mK)], Lf is the latent heat of fusion (J/Kg), rho_ice is the ice density, θ_ice is the soil ice content [m3/m3], t is the time, and z is the depth. The last term on the right side of Eq. (1) represents the latent heat of released/consumed during the phase change. The soil thermal conductivity follows the parameterization of Peters-Lidard (Ref) given by:

<img width="598" alt="eq2" src="https://user-images.githubusercontent.com/15165757/157315994-bca4d024-0288-4809-aeae-92e9012b4318.png">

where Ke , Ksat, and Kdry are the Kersten number, saturated, and dry soil thermal conductivities, respectively, and are defined in [Peters-Lidard et al., 1998](https://journals.ametsoc.org/view/journals/atsc/55/7/1520-0469_1998_055_1209_teostc_2.0.co_2.xml). The effective volumetric heat capacity, C, is calculated based on the respective fraction of each component (water, ice, air, and rock):

<img width="617" alt="eq3" src="https://user-images.githubusercontent.com/15165757/157316057-3aeffa5e-1b0e-4a29-b34a-12bda5d1e58e.png">

here, θ, θ_w, and θ_s are the total water content, the liquid water content and the maximum soil moisture content (porosity), respectively. The parameters  Cw, Cice, Cair, and Csoil are the volumetric heat capacities of water, ice, air, and soil, respectively.
The diffusion equation (Eq. 1) is discretized using the a fully implicit finite difference scheme and thus unconditionally stable.

Freezing-point depression equation:
Frozen soils not only impact the soil thermal and hydrological properties, but its presence strongly affects the partition of precipitation into surface runoff and infiltration leading to a strongly coupled surface and subsurface system. Due to the formation of impermeable soil over the winter period, infiltration of water significantly reduces during the spring freshet generating high peaks in the surface discharge.

In unfrozen conditions, the liquid water content is equal to the total water content in the soil. However, when soil freezes the liquid water content also depends on the soil temperature. The soil matric potential, Ψ_s [m], under the assumption that soil water potential remains in equilibrium with vapor pressure over pure ice when ice is present in the soil, and neglecting soil water osmotic potential, is given by:

<img width="607" alt="eq4" src="https://user-images.githubusercontent.com/15165757/157316409-3c51afda-e64b-4f6b-bf9e-efa50466cae4.png">	

where T and To  are the soil temperature and freezing temperatures, respectively (Ref), g is the acceleration of gravity [m/s^2] and Lf is the latent heat of fusion [J/Kg]. The soil matric potential, psi [m], as a function of soil moisture content is given by (ref: Clapp-Hornberger):

<img width="513" alt="eq5" src="https://user-images.githubusercontent.com/15165757/234057006-7b73ee81-4fb8-4320-bb55-91e7cf6f0b98.png">

where Ψ_s, θ_s, and θ_l are saturated soil matric potential [m], saturated soil moisture [-] (porosity), and soil liquid moisture [-], respectively. The exponent b is the Clapp-Hornberger parameter. According to the ‘freezing equals drying’ approximation (Ref), we equate Eq. (4) and Eq. (5) to obtain an expression for liquid water content:

<img width="573" alt="eq6" src="https://user-images.githubusercontent.com/15165757/234313360-8e07354e-877c-4c37-96bf-2dafe28ca0cf.png">

Equation (6) is called the “freezing-point depression equation” (Refs) and gives the maximum amount of liquid water (unfrozen soil moisture content) that can exist below the subfreezing temperature. The frozen soil moisture content (i, ice fraction) is given by:
θ_ice = θ_t- θ_l , where t is the total water content given by soil water retention curve: 

<img width="427" alt="eq7" src="https://user-images.githubusercontent.com/15165757/234313430-e9b64235-fe4d-47eb-9766-417957ffbb52.png">

