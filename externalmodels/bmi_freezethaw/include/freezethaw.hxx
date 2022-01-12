#ifndef FSH_INCLUDED
#define FSH_INCLUDED


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

using namespace std;


class Properties;
namespace freezethaw {
  
  class FreezeThaw{
  private:
    std::string config_file;
    std::string forcing_file;
    void InitializeArrays(void);
    
  public:
    int shape[3];
    double spacing[8];
    double origin[3];
    double time;
    double endtime;
    int nsteps;
    double dt;
    int nz;
    double lhf; // latent heat of fusion
    double tbot,ttop;
    double *Z; // depth
    double *Dz; // layer thickness
    double *ST; // soil temperature
    double *HC; //heat capacity
    double *TC; // thermal conductivity
    double *SMCT; // total soil moisture content
    double *SMCLiq; // liquid moisture content
    double *SMCIce; // ice moisture content
    double *GT; // ground/air temperature
    double *Time_;
    int opt_botb; //bottom boundary condition. 1 = zero flux, 2 = prescribed temperature
    int opt_topb; //top surface boundary condition. 1 = prescribed flux, 2 = prescribed temperature
    double smcmax; //porosity
    double bexp;  // pore size distribution [-], beta exponent on Clapp-Hornberger (1978)
    double psisat; // saturated capillary head (saturated moisture potential) [m]
    double ice_fraction;
    std::string ice_fraction_scheme;
    bool is_SMC_BMI_set;
    
    FreezeThaw();
    FreezeThaw(std::string config_file);
    
    void Advance();
    void SolveDiffusionEq();
    double GroundHeatFlux(double surfT);
    bool SolverTDMA(const vector<double> &a, const vector<double> &b, const vector<double> &c, const vector<double> &d, vector<double> &X);
    void PhaseChange();
    void ThermalConductivity();
    void SoilHeatCapacity();
    void SetLayerThickness();
    void InitFromConfigFile();
    double GetDt();
    
    std::vector<double> ReadVectorData(std::string key);
    void ReadForcingData(std::string key);
    void SetSMCBulk(); // set bulk moisture content per soil column
    ~FreezeThaw();
  };

  class Properties {
  private:
  public:
     //Volumetic heat capacity [j/m3/k]
    const double hcwater_; // water heat capacity
    const double hcice_; // ice heat capacity
    const double hcair_;  // air heat capacity
    const double hcsoil_; // rock/soil heat capacity
    const double tcice_;  // thermal conductiviyt of ice
    const double lhf_; // latent heat of fusion (j/kg)
    const double psisat_;  //Saturated matrix potential for soil type = silt loam
    const double grav_;
    const double tcwater_;// TC of water
    const double tcquartz_; // TC of Quartz
    const double quartz_; //loamy sand
    const double tcmineral_; // TC of other mineral 
    const double tfrez_;    // freezing/melting point (k)
    //const double smcmax_; // porosity (maximum soil moisture)
    const double bexp_; // Clap-Honnberger parameter
    const double wdensity_; // [kg/m3]
    Properties();
    ~Properties(){}
  };
};

#endif
