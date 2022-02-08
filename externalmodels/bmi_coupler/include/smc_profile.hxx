#ifndef SMCP_INCLUDED
#define SMCP_INCLUDED


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
    double water_table_m;
    double *water_table_prev_m;
    double *SMCT; // total soil moisture content


    double smcmax; //porosity
    double bexp;  // pore size distribution [-], beta exponent on Clapp-Hornberger (1978)
    double satpsi; // saturated capillary head (saturated moisture potential) [m]
    int nz;
    double D; //depth of the column/domain
    //    double *Z; // depth
    std::vector<double> Z;
    double *Dz; // layer thickness

    
    SMCProfile();
    SMCProfile(std::string config_file);
    
    void SetLayerThickness();
    void InitFromConfigFile();
    
    
    std::vector<double> ReadVectorData(std::string key);
    void ReadForcingData(std::string key);

    void SoilMoistureVerticalProfile();
      
    static const double grav;
    static const double wden;
    ~SMCProfile();
    
  };

};

#endif
