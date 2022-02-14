#ifndef SMCP_C_INCLUDED
#define SMCP_C_INCLUDED

#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "../include/smc_profile.hxx"
#define OK (1)

//const double smc_profile::SMCProfile::grav = 9.86;
//const double smc_profile::SMCProfile::wden = 1000.;

smc_profile::SMCProfile::
SMCProfile()
{
  this->shape[0] = 1;
  this->shape[1] = 1;
  this->shape[2] = 1;
  this->spacing[0] = 1.;
  this->spacing[1] = 1.;
  this->origin[0] = 0.;
  this->origin[1] = 0.;
  this->D =0.0;
  this->config_file = "";
  this->smcp_option= "";
  this->nz=0;
}

smc_profile::SMCProfile::
SMCProfile(std::string config_file)
{
  this->config_file = config_file;
  this->InitFromConfigFile();
  this->shape[0] = this->nz;
  this->shape[1] = 1;
  this->shape[2] = 1;
  this->spacing[0] = 1.;
  this->spacing[1] = 1.;
  this->origin[0] = 0.;
  this->origin[1] = 0.;
  this->InitializeArrays();
  
  SetLayerThickness(); // get soil layer thickness
  //  this->smcp_option= "";
}

void smc_profile::SMCProfile::
InitializeArrays(void)
{
  this->Dz = new double[nz];
  this->SMCT = new double[nz];
  this->SMCL = new double[nz];
  this->storage_m = new double[1];
  this->storage_change_m = new double[1];
  *this->storage_m = 0.0;
  *this->storage_change_m = 0.0;
  //  for (int i=0;i<nz;i++)
  //  this->SMCT[i] = 0.0;
}

void smc_profile::SMCProfile::
InitFromConfigFile()
{ 
  std::ifstream fp;
  fp.open(config_file);
  
  bool is_Z_set = false;
  bool is_Z_layers_set = false;
  bool is_smcmax_set = false;
  bool is_bexp_set = false;
  bool is_satpsi_set = false;
  bool is_wt_set = false;
  bool is_smcp_option_set = false;
  bool is_calc_linear_set = false;
  
  while (fp) {

    std::string key;
    std::getline(fp, key);
    
    int loc = key.find("=");
    std::string key_sub = key.substr(0,loc);
    
    if (key_sub == "Z") {
      std::string tmp_key = key.substr(loc+1,key.length());
      std::vector<double> vec = ReadVectorData(tmp_key);
      
      this->Z = new double[vec.size()];
      
      for (unsigned int i=0; i < vec.size(); i++)
	this->Z[i] = vec[i];
      
      this->nz = vec.size();
      this->D = this->Z[this->nz-1];
      is_Z_set = true;
      continue;
    }
    else if (key_sub == "Z_layers") {
      std::string tmp_key = key.substr(loc+1,key.length());
      std::vector<double> vec = ReadVectorData(tmp_key);
      
      this->Z_layers = new double[vec.size()];
      
      for (unsigned int i=0; i < vec.size(); i++)
	this->Z_layers[i] = vec[i];
      
      this->nz_layers = vec.size();
      this->D_layers = this->Z_layers[this->nz_layers-1];
      is_Z_layers_set = true;
      continue;
    }
    else if (key_sub == "soil_params.smcmax") {
      this->smcmax = std::stod(key.substr(loc+1,key.length()));
      is_smcmax_set = true;
      continue;
    }
    else if (key_sub == "soil_params.b") {
      this->bexp = std::stod(key.substr(loc+1,key.length()));
      assert (this->bexp > 0);
      is_bexp_set = true;
      continue;
    }
    else if (key_sub == "soil_params.satpsi") {  //Soil saturated matrix potential
      this->satpsi = std::stod(key.substr(loc+1,key.length()));
      is_satpsi_set = true;
      continue;
    }

    else if (key_sub == "soil_params.water_table") {
      double wt = std::stod(key.substr(loc+1,key.length()));
      this->water_table_m = new double[1];
      this->water_table_m[0] = this->D - wt;
      is_wt_set = true;
      continue;
    }

    else if (key_sub == "smcp_option") {  //Soil moisture profile option
      this->smcp_option = key.substr(loc+1,key.length());
      is_smcp_option_set = true;
      continue;
    }

    else if (key_sub == "calculated_linear") {  //Soil moisture profile option
      this->calc_linear = stod(key.substr(loc+1,key.length()));
      is_calc_linear_set = true;
      continue;
    }
  }
  fp.close();
  
  if (!is_Z_set) {
    std::stringstream errMsg;
    errMsg << "Z not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }

  if (!is_Z_layers_set) {
    std::stringstream errMsg;
    errMsg << "Z layers not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }
  
  if (!is_smcmax_set) {
    std::stringstream errMsg;
    errMsg << "smcmax not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }
  
  if (!is_bexp_set) {
    std::stringstream errMsg;
    errMsg << "bexp (Clapp-Hornberger's parameter) not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }
  
  if (!is_satpsi_set) {
    std::stringstream errMsg;
    errMsg << "satpsi not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }
  
  if (!is_wt_set) {
    std::cout<<"Warning! Water table location not provided, defualt is 1.9 m deep. \n";
    //initial water table location 1.9 m deep, if not provided in the config file
    this->water_table_m = new double[1];
    this->water_table_m[0] = this->D - 1.9; 
  }
  
  if (!is_smcp_option_set) {
    std::stringstream errMsg;
    errMsg << "soil moisture profile option not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }

  /*
  if(is_smcp_option_set) {
    if (this->calc_linear)
    //if (this->smcp_option != "calculated_linear")
    
  }*/
  
  // check if the size of the input data is consistent
  assert (this->nz >0);
  
}


std::vector<double> smc_profile::SMCProfile::
ReadVectorData(std::string key)
{
  int pos =0;
  std::string delimiter = ",";
  std::vector<double> value(0.0);
  std::string z1 = key;

  while (z1.find(delimiter) != std::string::npos) {
    pos = z1.find(delimiter);
    std::string z_v = z1.substr(0, pos);

    value.push_back(stod(z_v.c_str()));

    z1.erase(0, pos + delimiter.length());
    if (z1.find(delimiter) == std::string::npos)
      value.push_back(stod(z1));
  }

  return value;
}


void smc_profile::SMCProfile::
SetLayerThickness() {
  
  const int n_z = this->shape[0];

  Dz[0] = Z[0];
  for (int i=0; i<n_z-1;i++) {
    Dz[i+1] = Z[i+1] - Z[i];
    }
}

void smc_profile::SMCProfile::
SMPVertical()
{
  std::cout<<"SMP option "<<this->smcp_option<<"\n";
  if (this->smcp_option == "conceptual" || this->smcp_option == "Conceptual") {
    SMPFromConceptualReservoir();
  }
  else if (this->smcp_option == "calculated" || this->smcp_option == "Calculated") {
    SMPFromCalculatedReservoir();
  }
  else {
    std::stringstream errMsg;
    errMsg << "Soil moisture profile OPTION provided in the config file "<< config_file << " is " << this->smcp_option<< ", which should be either \'concepttual\' or \'calculated\' " <<"\n";
    throw std::runtime_error(errMsg.str());

  }
}

// given bulk soil moisture quantity, distribute vertically using CFE conceptual storage
void smc_profile::SMCProfile::
SMPFromConceptualReservoir()
{
  double lam=1.0/this->bexp; // pore distribution index
  double hb=this->satpsi * 100.; //7.82;  //[cm] //soil_res->hp;     
  double D= this->D * 100.;
  double z1= *this->water_table_m * 100; //previous water table location, [cm]
  
  double Vmax=D* this->smcmax;
  double beta=1.0-lam;
  double alpha=pow(hb,lam)/beta;
  double tol=0.000001;
  double phi = this->smcmax;
  double V1 = 100.0 * (*this->storage_m - *this->storage_change_m);  //change in soil moisture

  double V2= 100.0 * (*this->storage_m);  /* start-up condition before adding any water */
  int count = 0;
  
  if(V2>=Vmax) {
    for(int j=0;j<this->nz;j++)
      this->SMCT[j] = phi;
    return;
  }
   
   double diff=1000.0;

   double f, z2, z2new, df_dz2;

   z2=z1;  /* first guess, use Newton-Raphson to find new Z2 */
   do {
     count++;
    
     if(count>15000) {
       throw std::runtime_error("No convergence loop count: after 15000 iterations!");
     }
     
     f=phi*(z2-z1) + alpha * phi * (std::pow((D-z2),beta)-std::pow((D-z1),beta)) - (V2-V1);

     df_dz2=phi - alpha*phi*beta*std::pow((D-z2),(beta-1.0));
     
     z2new=z2-f/df_dz2*1.0;
     
     diff=z2new-z2;

     z2=z2new;     
     
   } while (std::fabs(diff)>tol);
   
   z1=z2;  // reset to new water table elevation value
   *this->water_table_m = z1/100.;

   /* get a high resolution curve */
   int z_hres = 1001;
   double *smct_temp = new double[z_hres];
   double *z_temp = new double[z_hres];
   double dz1 = hb;
   double dz2 = this->Z[this->nz-1]*100.0 - z1;
   
   for (int i=0;i<z_hres;i++) {
     smct_temp[i] = std::pow((hb/dz1),lam)*phi;
     z_temp[i] = z1  + dz1;
     dz1 += dz2/(z_hres-1);
   }
   
   // mapping the updated soil moisture curve to the heat conduction discretization depth (Dz)
   for (int i=0; i<this->nz; i++) {
     for (int j=0; j<z_hres; j++) {
       if (z_temp[j]  >= (D - this->Z[i]*100) ) {
	 this->SMCT[i] = smct_temp[j];
	 break;
       }
     }
   }
   
}

// given layered-based soil moisture (calculated storage such as LGAR), distribute vertically to match soil freeze thaw resolution
void smc_profile::SMCProfile::
SMPFromCalculatedReservoir()
{

  double lam=1.0/this->bexp; // pore distribution index
  double hb=this->satpsi; //7.82;  //[cm] //soil_res->hp;
  double phi = this->smcmax;
  
  double D = this->D; //6.0; // depth of the SFT model
  double Dl = this->D_layers; //2.0; // depth of the layer model
  int nlayer = this->nz_layers;
  
  std::vector<double> z_layers_n(1,0.0);
  
  for (int i=0; i <this->nz_layers; i++)
    z_layers_n.push_back(Z_layers[i]);

  
  std::vector<double> z_res; // from SFT
  int sft_ncells = 60;
  
  for (int i=0; i < sft_ncells; i++)
    z_res.push_back((i+1)*0.1);

  std::vector<double> smc_sft;
  int c =0;
  double delta =0.0;
  
  // piece-wise constant (vertically)
  if (!this->calc_linear) {
    bool layers_flag=true;
    for (int i=0; i < sft_ncells; i++){

      if (z_res[i] < Z_layers[c]) {
	smc_sft.push_back(SMCL[c]);
      }
      else if (z_res[i] < Dl) {
	double v_avg = 0.5*(SMCL[c] + SMCL[c+1]); // interface of layers
	smc_sft.push_back(v_avg);
	c++;
      }
      else {
	double zz = D - z_res[i];
	double theta = delta + std::pow((hb/zz),lam)*phi;
      
	if (layers_flag) { // to make sure the moisture curve below Dl start at moisture value at Dl for continuity
	  delta = 0.0 - theta;
	  theta = theta + delta;
	  layers_flag=false;
	}
	
	double theta1 = SMCL[nlayer-1] +  theta; //pow((hb/z2),lam)*phi
	double theta2 = std::min(phi, theta1);
	if (z_res[i] > hb)
	  smc_sft.push_back(theta2);
	else
	  smc_sft.push_back(phi);
      }
      std::cout<<"cons. "<<smc_sft[i]<<","<<" \n";

    }
  }
  else {
    // linearly interpoloted
    bool layers_flag=true;
    double t_v=0.0;
    
    smc_sft.push_back(SMCL[0]); // SFT first cell gets the top-layer soil moisture
    
    for (int i=1; i < sft_ncells; i++){
      
      if (z_res[i] <= Z_layers[c]) {
	
	if (c == nlayer-1)
	  t_v = LinearInterpolation(z_layers_n[c], z_layers_n[c+1], SMCL[c], SMCL[c], z_res[i]);// interface of layers
	else
	  t_v = LinearInterpolation(z_layers_n[c], z_layers_n[c+1], SMCL[c], SMCL[c+1], z_res[i]);
	smc_sft.push_back(t_v);
	
      }
      else if (z_res[i] < Dl) {
	
	c++;
	
	if (c == nlayer-1)
	  t_v = LinearInterpolation(z_layers_n[c], z_layers_n[c+1], SMCL[c], SMCL[c], z_res[i]);// interface of layers
	else
	  t_v = LinearInterpolation(z_layers_n[c], z_layers_n[c+1], SMCL[c], SMCL[c+1], z_res[i]);// interface of layers
	smc_sft.push_back(t_v);
	
      }
      else {
	
	double zz = D - z_res[i];
	
	double theta = delta + std::pow((hb/zz),lam)*phi;
	
	if (layers_flag) { // to make sure the moisture curve below Dl start at moisture value at Dl for continuity
	  delta = 0.0 - theta;
	  theta = theta + delta;
	  layers_flag = false;
	  
	}
	double theta1 = SMCL[nlayer-1] +  theta;
	double theta2 = std::min(phi, theta1);
	if (z_res[i] > hb)
	  smc_sft.push_back(theta2);
	else
	  smc_sft.push_back(phi);
	
      }
       std::cout<<smc_sft[i-1]<<",";
      
      }
      
    }
}

double smc_profile::SMCProfile::
LinearInterpolation(double z1, double z2, double t1, double t2, double z)
{

  double m = (t2 - t1) / (z2 - z1);
  return t1 + m * (z - z1);

}

smc_profile::SMCProfile::
~SMCProfile()
{}

#endif
