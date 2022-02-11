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

    double smcmax; //porosity
    double bexp;  // pore size distribution [-], beta exponent on Clapp-Hornberger (1978)
    double satpsi; // saturated capillary head (saturated moisture potential) [m]
    int nz;
    double D; //depth of the column/domain
    //std::vector<double> Z;
    double *Z;
    double *Dz; // layer thickness
    std::string smcp_option;
    
    SMCProfile();
    SMCProfile(std::string config_file);
    
    void SetLayerThickness();
    void InitFromConfigFile();
    
    std::vector<double> ReadVectorData(std::string key);
    void ReadForcingData(std::string key);

    void SMPVertical(); //SMP -> Soil Moisture Profile

    void SMPFromConceptualReservoir();

    void SMPFromCalculatedReservoir();
      
    //    static const double grav;
    //    static const double wden;
    ~SMCProfile();
    
  };

};

#endif
