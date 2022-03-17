## Building the code to run/test examples within the pseudo-framework (not the ngen framework!!)
### Steps to build soil freeze-thaw (SFT) model coupled with cfe and soil moisture profiles
### Setting up (cloning) external repos
- git clone https://github.com/NOAA-OWP/SoilFreezeThaw && cd SoilFreezeThaw (**for now: git checkout ajk/sft_only**)
- git clone https://github.com/NOAA-OWP/cfe 
- git checkout [cfe_soilfreezethaw](https://github.com/NOAA-OWP/cfe/tree/cfe_soilfreezethaw) (note: coupling SFT with CFE using the pseudo-framework requires this special cfe branch)
- git clone https://github.com/NOAA-OWP/SoilMoistureProfiles smc_coupler
- mkdir build && cd build
- cmake -DCMAKE_INSTALL_PREFIX=\`pwd\` -DCMAKE_BUILD_TYPE=Debug ../
- make && cd ..
- [run_framework.sh](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/sft_only/run_framework.sh)



# Introduction of Soil Freeze-thaw model

The diffusion equation is used to simulate the transport of energy in the soil.

<img width="615" alt="eq1" src="https://user-images.githubusercontent.com/15165757/157314534-0ef6e5ea-4dad-4be5-aaac-888ca139cbc7.png">

where T is the soil temperature, C is the volumetric heat capacity [J/(m3K)], θ is the volumetric soil moisture content [m3/m3], K is the bulk thermal conductivity [W/(mK)], Lf is the latent heat of fusion (J/Kg), ice is the ice density, ice is the soil ice content [m3/m3], t is the time, and z is the depth. The last term on the right side of Eq. (1) represents the latent heat of released/consumed during the phase change. The soil thermal conductivity follows the parameterization of Peters-Lidars (Ref) given by:

<img width="598" alt="eq2" src="https://user-images.githubusercontent.com/15165757/157315994-bca4d024-0288-4809-aeae-92e9012b4318.png">

where Ke , Ksat, and Kdry are the Kersten number, dry, and saturated soil thermal conductivities, respectively, and are defined in (Ref). The effective volumetric heat capacity, C, is calculated based on the respective fraction of each component (water, ice, air, and rock):

<img width="617" alt="eq3" src="https://user-images.githubusercontent.com/15165757/157316057-3aeffa5e-1b0e-4a29-b34a-12bda5d1e58e.png">

here, θ, θ_w, and s are the total water content, the liquid water content and the maximum soil moisture content (porosity), respectively. The parameters  Cw, Cice, Cair, and Csoil are the volumetric heat capacities of water, ice, air, and soil, respectively. Table?
The diffusion equation (Eq. 1) is discretized using the Crank-Nicolson scheme (Ref), which is an implicit finite difference scheme and thus unconditionally stable.

Freezing-point depression equation 
Frozen soils not only impact the soil thermal and hydrological properties, but its presence strongly affects the partition of precipitation into surface runoff and infiltration leading to a strongly coupled surface and subsurface system. Due to the formation of impermeable soil over the winter period, infiltration of water significantly reduces during the spring freshet generating high peaks in the surface discharge.

In unfrozen conditions, the liquid water content is equal to the total water content in the soil. However, when soil freezes the liquid water content also depends on the soil temperature. The soil matric potential,  [m], under the assumption that soil water potential remains in equilibrium with vapor pressure over pure ice when ice is present in the soil, and neglecting soil water osmotic potential, is given by:

<img width="607" alt="eq4" src="https://user-images.githubusercontent.com/15165757/157316409-3c51afda-e64b-4f6b-bf9e-efa50466cae4.png">	

where T and T0  are the soil temperature and freezing temperatures, respectively (Ref), g is the acceleration of gravity [m s-2] and Lf is the latent heat of fusion [J Kg-1]. The soil matric potential,  [m], as a function of soil moisture content is given by (ref: Clapp-Hornberger):

<img width="594" alt="eq5" src="https://user-images.githubusercontent.com/15165757/157317205-a2b1e455-8714-4113-8c56-f7a2ef55533c.png">

where s, s, and l are saturated soil matric potential [m], saturated soil moisture [-] (porosity), and soil liquid moisture [-], respectively. The exponent b is the Clapp-Hornberger parameter. According to the ‘freezing equals drying’ approximation (Ref), we equate Eq. (4) and Eq. (5) to obtain an expression for liquid water content:

<img width="626" alt="eq6" src="https://user-images.githubusercontent.com/15165757/157317016-3b8c0ad1-1f4b-4370-9773-c97aba0d0ee5.png">

Equation (6) is called the “freezing-point depression equation” (Refs) and gives the maximum amount of liquid water (unfrozen soil moisture content) that can exist below the subfreezing temperature. The frozen soil moisture content (i, ice fraction) is given by:
i= t- θl , where t is the total water content given by soil water retention curve: 
t= θss-b for T>T0.

