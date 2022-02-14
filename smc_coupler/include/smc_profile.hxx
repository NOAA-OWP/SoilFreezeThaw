#ifndef SMCP_H_INCLUDED
#define SMCP_H_INCLUDED


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

using namespace std;

namespace smc_profile {
  
  class SMCProfile{
  private:
    std::string config_file;
    void InitializeArrays(void);
    
  public:
    int shape[3];
    double spacing[8];
    double origin[3];

    double *storage_m;
    double *storage_change_m;
    double *water_table_m;
    double *SMCT; // total soil moisture content
    double *SMCL; // total layered-soil moisture content
    
    double smcmax; //porosity
    double bexp;  // pore size distribution [-], beta exponent on Clapp-Hornberger (1978)
    double satpsi; // saturated capillary head (saturated moisture potential) [m]
    int nz;
    int nz_layers;
    double D; //depth of the column/domain
    double D_layers; // depth of the last layer
    //std::vector<double> Z;
    double *Z;
    double *Z_layers;
    double *Dz; // layer thickness
    std::string smc_profile;
    std::string smc_profile_option;
    int *smc_profile_option_bmi;
    
    SMCProfile();
    SMCProfile(std::string config_file);
    
    void SetLayerThickness();
    void InitFromConfigFile();
    
    std::vector<double> ReadVectorData(std::string key);
    void ReadForcingData(std::string key);

    void SMPVertical(); //SMP -> Soil Moisture Profile

    void SMPFromConceptualReservoir();

    void SMPFromLayeredReservoir();

    double LinearInterpolation(double z1, double z2, double t1, double t2, double z);
    //    static const double grav;
    //    static const double wden;
    ~SMCProfile();
    
  };

};

#endif
