/*
  @author Ahmad Jan (ahmad.jan@noaa.gov)
  @date   September 2021
*/

/*
  The code simulates the transport of energy in the soil using the 1D diffusion equation (for more details see README on github repo (https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/README.md).
  Crank-Nicolson scheme is used to discretized the equation
  OUTPUTS : soil_temperature, ice_fraction_schaake, ice_fraction_xinan
  INPUTS  : Rest of the parameters are inputs either through a BMI or config file
  Soil thermal conductivity and volumetric heat capacity use empirical models
  The model is coupled to the surface through ground surface temperature
  The model is also coupled to soil moisture profile module to provide dynamic vertical distribution of soil moisture

  @param soil_depth [m]                    : depth of the computational domain 
  @param latent_heat_fusion         [J/kg] : latent heat of fusion
  @param ground_temp_const          [K]    : constant ground temperature (an option for the surface boundary condition)
  @param end_time                   [s]    : end time of the simulation, input options [second, hour, day]
  @param dt                         [s]    : timestep, input options [second, hour, day]
  @param ncells                     [-]    : number of cells in the discretized soil column
  @param bottom_boundary_temp_const [K]    : temperature at the bottom boundary of the domain
  @param topboundary_temp_const     [K]    : temperature at the top boundary of the domain (ground temperature)
  @param soil_z                     [m]    : soil discretization, depth from the surface
  @param soil_dz                    [m]    : soil discretization thickness, thickness of cells
  @param soil_temperature           [K]    : soil temperature profile
  @param heat_capacity              [J/(m3 K)] : volumetric heat capacity (specific heat capacity * density)
  @param thermal_conductivity       [W/(mK)]   : soil bulk thermal conductivity
  @param soil_moisture_content      [-]    : total (ice+water) soil moisture content (1D profile)
  @param soil_liquid_content        [-]    : portion of liquid in the total soil moisture content (1D profile)
  @param soil_ice_content           [-]    : portion of ice in the total soil moisture content (1D profile)
  @param soil_ice_fraction          [-]    : fraction of soil moisture that is ice (scalar)
  @param option_bottom_boundry      [-]    : option for bottom boundary condition. 1 = zero geothermal flux, 2 = prescribed temperature
  @param option_top_boundary        [-]    : top surface boundary condition. 1 = prescribed flux, 2 = prescribed temperature
  @param smcmax                     [-]    : maximum soil moisture content (porosity)
  @param bb                         [-]    : pore size distribution, beta exponent in Clapp-Hornberger (1978) function
  @param satpsi                     [m]    : saturated capillary head (saturated moisture potential, capillary fringe thickness)
  @param ice_fraction_schaake       [-]    : ice fraction based on Schaake runoff scheme (computes volume of frozen water)
  @param ice_fraction_xinan         [-]    : ice fraction based on Xinanjiang runoff scheme (based on ice content in the top cell)
  @param ice_fraction_scheme        [-]    : option set from the config file for the runoff scheme (Schaake=1, Xinanjiang=2)
  @param ice_fraction_scheme_bmi    [-]    : option set through a bmi
  @param is_soil_moisture_bmi_set   [-]    : if not standalone, soil moisture is set through SoilMoistureProfiles bmi
  @param quartz                     [-]    : quartz content used in the thermal conductivity model
*/

#ifndef SFT_H_INCLUDED
#define SFT_H_INCLUDED


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

using namespace std;

extern string verbosity;

class Properties;

namespace soilfreezethaw {
  
  class SoilFreezeThaw {
  private:
    std::string config_file;
    void InitializeArrays(void);
    
  public:
    int    shape[3];
    double spacing[8];
    double origin[3];
    double time;
    double endtime;
    double dt;
    int ncells;
    double soil_depth; 
    double latent_heat_fusion;
    double bottom_boundary_temp_const;
    double top_boundary_temp_const;
    double ground_temp_const;
    double *soil_z = NULL;
    double *soil_dz = NULL;
    double *soil_temperature = NULL;
    double *heat_capacity = NULL;
    double *thermal_conductivity = NULL;
    double *soil_moisture_content = NULL;
    double *soil_liquid_content = NULL;
    double *soil_ice_content = NULL;
    double soil_ice_fraction;
    double ground_temp;
    int    option_bottom_boundary;
    int    option_top_boundary; 
    double smcmax; 
    double bb;
    double satpsi;
    double quartz;
    double ice_fraction_schaake;
    double ice_fraction_xinan;
    int    ice_fraction_scheme_bmi;
    bool   is_soil_moisture_bmi_set;
    
    
    std::string ice_fraction_scheme;
    enum SurfaceRunoffScheme{Schaake=1, Xinanjiang=2}; // surface runoff schemes

    /*
      input var changes when running SFT in coupled mode and in standalone mode in the ngen framework. When standalone, soil_moisture_profile is not needed, so we have to setup input var names dynamically; only needed for ngen framework
     */
    std::vector<std::string>* input_var_names_model; 
    
    SoilFreezeThaw();
    SoilFreezeThaw(std::string config_file);
    
    void Advance();
    void SolveDiffusionEquation();
    double GroundHeatFlux(double surfT);

    /* Tridiagonal matrix solver */
    bool SolverTDMA(const vector<double> &a, const vector<double> &b, const vector<double> &c, const vector<double> &d, vector<double> &X); 

    /* Phase change module using freezing-point depression model */
    void PhaseChange();

    /* Thermal conductivity module using Peters-Lidard scheme*/
    void ThermalConductivity();

    /* computes volumetic heat capacity*/
    void SoilHeatCapacity();

    /* returns cells thickness*/
    void SoilCellsThickness();
    
    void InitFromConfigFile(std::string config_file);
    double GetDt();
    
    std::vector<double> ReadVectorData(std::string key);

    /* computes surface runoff-based ice fraction*/
    void ComputeIceFraction();

    // method retuns dynamically allocated input variable names
    std::vector<std::string>* InputVarNamesModel();
    
    ~SoilFreezeThaw();
  };

  // class to contain constant variables
  class Properties {
  private:
  public:
    const double hcwater_;  // water heat capacity
    const double hcice_;    // ice heat capacity
    const double hcair_;    // air heat capacity
    const double hcsoil_;   // rock/soil heat capacity
    const double grav_;     // gravity
    const double tfrez_;    // freezing/melting point (k)
    const double wdensity_; // liquid density [kg/m3]
    Properties();
    ~Properties(){}
  };
};

#endif
