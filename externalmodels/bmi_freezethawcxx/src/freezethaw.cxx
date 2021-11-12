#ifndef FSC_INCLUDED
#define FSC_INCLUDED

#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "../include/freezethaw.hxx"
#define OK (1)


freezethaw::FreezeThaw::
FreezeThaw()
{
  this->endtime = 10.;
  this->time = 0.;
  this->shape[0] = 1;
  this->shape[1] = 1;
  this->shape[2] = 1;
  this->spacing[0] = 1.;
  this->spacing[1] = 1.;
  this->origin[0] = 0.;
  this->origin[1] = 0.;
  this->dt = 3600;
  this->lhf = 0.3336E06;
  this->ttop = 260.;
  this->tbot = 272.;
  this->opt_botb = 1;
  this->opt_topb = 2;
  this->smcliq_bulk = 0.;
  this->smcice_bulk = 0.;
  this->smct_bulk = 0.;
  this->forcing_file= " ";
  this->nsteps=0;
}

freezethaw::FreezeThaw::
FreezeThaw(std::string config_file)
{

  this->config_file = config_file;
  this->forcing_file= " ";
  this->lhf = 0.3336E06;
  this->ttop = 260.;
  this->tbot = 272.;
  this->opt_botb = 1;
  this->opt_topb = 2; // 1: constant temp, 2: from a file

  this->InitFromConfigFile();

  this->shape[0] = this->nz;
  this->shape[1] = 1;
  this->shape[2] = 1;
  this->spacing[0] = 1.;
  this->spacing[1] = 1.;
  this->origin[0] = 0.;
  this->origin[1] = 0.;

  this->smcliq_bulk = 0.;
  this->smcice_bulk = 0.;
  this->smct_bulk = 0.;
  this->InitializeArrays();
  SetLayerThickness(); // get soil layer thickness
  SetSMCBulk();
  this->time = 0.;
  this->nsteps = 0;
}


void freezethaw::FreezeThaw::
InitializeArrays(void)
{
  this->TC = new double[nz];
  this->HC = new double[nz];
  this->Dz = new double[nz];
  this->SMCIce = new double[nz];
  
  for (int i=0;i<nz;i++)
    this->SMCIce[i] = this->SMCT[i] - this->SMCLiq[i];
}

int freezethaw::FreezeThaw::
InitFromConfigFile()
{ 
  std::ifstream fp;
  fp.open(config_file);
  
  while (fp) {
    
    std::string key;
    std::getline(fp, key);

    int loc = key.find("=");
    std::string key_sub = key.substr(0,key.find("="));

    if (key_sub == "forcing_file") {
      //this->forcing_file = key.substr(loc+1,key.length());
      std::string tmp_key = key.substr(loc+1,key.length());
      this->ReadForcingData(tmp_key);
      continue;
    }
    if (key_sub == "end_time_d") {
      this->endtime = std::stod(key.substr(loc+1,key.length()));
      this->endtime *= 86400;
      continue;
    }
    if (key_sub == "dt_s") {
      this->dt = std::stod(key.substr(loc+1,key.length()));
      continue;
    }
    if (key_sub == "Z") {
      std::string tmp_key = key.substr(loc+1,key.length());
      std::vector<double> vec = ReadVectorData(tmp_key);
      this->Z = new double[vec.size()];
      for (int i=0; i < vec.size(); i++)
	this->Z[i] = vec[i];
      this->nz = vec.size();
      continue;
    }
    if (key_sub == "nz") {
      int nz_t = std::stod(key.substr(loc+1,key.length()));
      //      assert (nz_t == this->nz);
      continue;
    }
    if (key_sub == "soil_params.smcmax") {
      this->smcmax = std::stod(key.substr(loc+1,key.length()));
      continue;
    }
    if (key_sub == "soil_temperature") {
      std::string tmp_key = key.substr(loc+1,key.length());
      std::vector<double> vec = ReadVectorData(tmp_key);
      this->ST = new double[vec.size()];
      for (int i=0; i < vec.size(); i++)
	this->ST[i] = vec[i];
      continue;

    }
    if (key_sub == "soil_total_moisture_content") {
      std::string tmp_key = key.substr(loc+1,key.length());
      std::vector<double> vec = ReadVectorData(tmp_key);
      this->SMCT = new double[vec.size()];
      for (int i=0; i < vec.size(); i++) {
	this->SMCT[i] = vec[i];
      }
      continue;
    }
    if (key_sub == "soil_liquid_moisture_content") {
      std::string tmp_key = key.substr(loc+1,key.length());
      std::vector<double> vec = ReadVectorData(tmp_key);
      this->SMCLiq = new double[vec.size()];
      for (int i=0; i < vec.size(); i++) {
	//	assert (this->SMCT[i] >= vec[i]);
	this->SMCLiq[i] = vec[i];
      }
      continue;
    }
  }

  fp.close();
  return 1;
}


std::vector<double> freezethaw::FreezeThaw::
ReadVectorData(std::string key)
{
  int pos =0;
  std::string delimiter = ",";
  std::vector<double> z_value(0.0);
  std::string z1 = key;
  
  while ((pos = z1.find(delimiter)) != std::string::npos) {
    //std::ostringstream z_vv;
    //z_vv << std::setprecision(8) << z1.substr(0, pos);
    std::string z_v = z1.substr(0, pos);
    z_value.push_back(stod(z_v.c_str()));

    z1.erase(0, pos + delimiter.length());
    if (z1.find(delimiter) == std::string::npos)
      z_value.push_back(stod(z1));
  }
  
  return z_value;
}

void freezethaw::FreezeThaw::
ReadForcingData(std::string forcing_file)
{
  std::ifstream fp;
  fp.open(forcing_file);
  std::vector<double> Time_v(0.0);
  std::vector<double> GT_v(0.0);
  std::vector<string> vars;
  std::string line, cell;

  //read first line of strings which contains forcing variables names.
  std::getline(fp, line);
  std::stringstream lineStream(line);
  while(std::getline(lineStream,cell, ',')) {
    vars.push_back(cell);    
  }

  int len_v = vars.size(); // number of forcing variables + time
  int count = 0;
  while (fp) {
    std::getline(fp, line);
    std::stringstream lineStream(line);
    while(std::getline(lineStream,cell, ',')) {
      
      if (count % len_v == 0) {
	Time_v.push_back(stod(cell));
	count +=1;
	continue;
      }
      if (count % len_v == 6) {
	// 6 is the air temperature column, needs to be fixed
	GT_v.push_back(stod(cell));
	count +=1;
	continue;
      }
      count +=1;
    }

  }

  int size_v = Time_v.size();

  this->Time_ = new double[size_v];
  this->GT = new double[size_v];

  
  for (int i=0; i<size_v; i++) {
    this->Time_[i] = Time_v[i];
    this->GT[i] = GT_v[i];
  }

}

void freezethaw::FreezeThaw::
SetSMCBulk()
{
  double val = 0;
  for (int i =0; i < nz; i++) val += this->SMCT[i];
  this->smct_bulk = val;

  val = 0;
  for (int i =0; i < nz; i++) val += this->SMCLiq[i];
  this->smcliq_bulk = val;

  val = 0;
  for (int i =0; i < nz; i++) val += this->SMCIce[i];
  this->smcice_bulk = val;

  this->ice_fraction = this->smcice_bulk / this->smct_bulk;
  assert (this->ice_fraction <= 1.0);
}
  
double freezethaw::FreezeThaw::
GetDt()
{
  return this->dt;
}

void freezethaw::FreezeThaw::
UpdateSMCDistribution()
{
  double d; // depth from the top surface to the last soil layer
  double d_wt_dry; // depth of the water table for drier conditions
  double d_wt_cur; // depth of the current water table
  double smc_top;  // smc in the soil layers (entire domain) [0,d]
  double smc_bot;  // smc between d and d_wt_cur depths
  double smc_max;  // maximum smc between d and d_wt_dry depths
  double mass_w;   // water mass added or removed
  double smc_total; // smc_top + smc_bot + mass_w
}


void freezethaw::FreezeThaw::
Advance()
{
  // Update Thermal conductivities, note the soil heat flux update happens in the PhaseChange module, so no need to update here
  ThermalConductivity(); // initialize thermal conductivities

  // call Diffusion Eq. first to get the updated soil temperatures
  SolveDiffusionEq();
  
  // call phase change to get water/ice partition for the updated temperature
  PhaseChange();

  this->time += this->dt;

  SetSMCBulk();
  this->nsteps += 1;
}

double freezethaw::FreezeThaw::
GroundHeatFlux(double surfT)
{  
  if (opt_topb == 1) {
    double ghf = - TC[0] * (surfT  - ttop) / Z[0];   //temperature specified as constant
    return ghf; 
  }
  else if (opt_topb == 2) {
    double ghf = - TC[0] * (surfT  - GT[this->nsteps]) / Z[0];  //temperature from a file
    return ghf; 
  }
  else
    return 0;
}

//*********************************************************************************
// Solve, using the Thomas Algorithm (TDMA), the tri-diagonal system              *
//     a_i X_i-1 + b_i X_i + c_i X_i+1 = d_i,     i = 0, n - 1                    *
//                                                                                *
// Effectively, this is the n x n matrix equation.                                *
// a[i], b[i], c[i] are the non-zero diagonals of the matrix and d[i] is the rhs. *
// a[0] and c[n-1] aren't used.                                                   *
//*********************************************************************************
bool freezethaw::FreezeThaw::
SolverTDMA(const vector<double> &a, const vector<double> &b, const vector<double> &c, const vector<double> &d, vector<double> &X ) {
   int n = d.size();
   vector<double> P( n, 0 );
   vector<double> Q( n, 0 );
   X = P;
   // Forward pass
   double denominator = b[0];

   P[0] = -c[0]/denominator;
   Q[0] =  d[0]/denominator;

   for (int i = 1; i < n; i++) {
     denominator = b[i] + a[i] * P[i-1];

     if ( std::abs(denominator) < 1e-20 ) return false;
     
     P[i] =  -c[i]/denominator;
     Q[i] = (d[i] - a[i] * Q[i-1])/denominator;
   }
   
   // Backward substiution
   X[n-1] = Q[n-1];
   for (int i = n - 2; i >= 0; i--)
     X[i] = P[i] * X[i+1] + Q[i];
   
   return true;
}
  
void freezethaw::FreezeThaw::
SolveDiffusionEq ()
{
    int i;
    const int nz = this->shape[0];

    // local variables
    std::vector<double> Flux(nz);
    std::vector<double> AI(nz);
    std::vector<double> BI(nz);
    std::vector<double> CI(nz);
    std::vector<double> RHS(nz);
    std::vector<double> Lambd(nz);
    std::vector<double> X(nz);
    double h, botflux;
    // compute matrix coefficient using Crank-Nicolson discretization scheme
    for (i=0;i<nz; i++) {
      if (i == 0) {
	h = Z[i];
	double dtdz  = (ST[i+1] - ST[i])/ h;
	Lambd[i] = dt/(4* h * lhf);
	double ghf = this->GroundHeatFlux(ST[0]);
	Flux[i] = Lambd[i] * (TC[i] * dtdz + ghf);
      }
      else if (i < nz -1) {
	h = Z[i+1] - Z[i];
        double h1 = Z[i] - Z[i-1];
	Lambd[i] = dt/(4*h* lhf);
	double a_ = - Lambd[i] * TC[i] /h1;
        double c_ = - Lambd[i] * TC[i+1] /h;
	double b_ = 1 + a_ + c_;
	Flux[i] = -a_ * ST[i-1] + b_ * ST[i] - c_ * ST[i+1];
      }
      else if (i == nz-1) {
	h = (Z[i] - Z[i-1]);
	Lambd[i] = dt/(4* h* lhf);
	if (opt_botb == 1) 
	  botflux = 0.;
	else if (opt_botb == 2) {
	  double dtdz1 = (ST[i] - tbot) / ( Z[i] - Z[i-1]);
	  botflux  = - TC[i] * dtdz1;
	}
	double dtdz = (ST[i] - ST[i-1] )/ h;
	Flux[i]  = Lambd[i] * (-TC[i]*dtdz  + botflux);
      }
    }
    // put coefficients in the corresponding vectors A,B,C, RHS
    for (i=0; i<nz;i++) {
      if (i == 0) {
	AI[i] = 0;
	CI[i] = - Lambd[i] *TC[i]/Z[i];
	BI[i] = 1 - CI[i];
      }
      else if (i < nz-1) {
	AI[i] = - Lambd[i] * TC[i]/(Z[i] - Z[i-1]);
	CI[i] = - Lambd[i] * TC[i+1]/(Z[i+1] - Z[i]);
	BI[i] = 1 - AI[i] - CI[i];
      }
      else if (i == nz-1) { 
	AI[i] = - Lambd[i] * TC[i]/(Z[i] - Z[i-1]);
	CI[i] = 0;
	BI[i] = 1 - AI[i];
      }
      RHS[i] = Flux[i];
    }

    // add the previous timestep ST to the RHS at the boundaries
    for (i=0; i<nz;i++) {
      if (i ==0)
	RHS[i] = ST[i] + RHS[i];
      else if (i == nz-1)
	RHS[i] = ST[i] + RHS[i];
      else
	RHS[i] = RHS[i];
    }

    SolverTDMA(AI, BI, CI, RHS, X);

    std::copy(X.begin(), X.end(), this->ST);
}


void freezethaw::FreezeThaw::
ThermalConductivity() {
  Properties prop;
  const int n_z = this->shape[0];
  
  for (int i=0; i<n_z;i++) {
    double sat_ratio = SMCT[i]/ prop.smcmax_;
    
    //TC of solids Eq. (10) Peters-Lidard
    double tc_solid = pow(prop.tcquartz_,prop.quartz_) * pow(prop.tcmineral_, (1. - prop.quartz_));

    //SATURATED THERMAL CONDUCTIVITY
    
    //UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
    double x_unfrozen = SMCLiq[i] / SMCT[i]; // (phi * Sliq) / (phi * sliq + phi * sice) = sliq/(sliq+sice) 
    
    double xu = x_unfrozen * prop.smcmax_; // unfrozen volume fraction
    double tc_sat = pow(tc_solid,(1. - prop.smcmax_)) * pow(prop.tcice_, (prop.smcmax_ - xu)) * pow(prop.tcwater_,xu);
    
    //DRY THERMAL CONDUCTIVITY
    double gammd = (1. - prop.smcmax_)*2700.; // dry density
    double tc_dry = (0.135* gammd+ 64.7)/ (2700. - 0.947* gammd);
    
    // Kersten Number
    // for frozen soil
    double KN;
    if (SMCLiq[i] + 0.001 > SMCT[i])
      KN = sat_ratio;
    else
      KN = sat_ratio > 0.1 ? log10(sat_ratio) + 1. : 0.0;
    
    // Thermal conductivity
    TC[i] = KN * (tc_sat - tc_dry) + tc_dry;
    
  }
}

void freezethaw::FreezeThaw::
SoilHeatCapacity() {
  Properties prop;
  const int n_z = this->shape[0];
  //def soil_heat_capacity(domain, prop,SMC, SOLIQ, HCPCT, SMCMAX):
  for (int i=0; i<n_z;i++) {
    double sice = SMCT[i] - SMCLiq[i];
    HC[i] = SMCLiq[i]*prop.hcwater_ + sice*prop.hcice_ + (1.0-prop.smcmax_)*prop.hcsoil_ + (prop.smcmax_-SMCT[i])*prop.hcair_;
  }

}

void freezethaw::FreezeThaw::
SetLayerThickness() {
  const int n_z = this->shape[0];

  Dz[0] = Z[0];
  for (int i=0; i<n_z-1;i++) {
    Dz[i+1] = Z[i+1] - Z[i];
  }
}

void freezethaw::FreezeThaw::
PhaseChange() {
  Properties prop;
  const int n_z = this->shape[0];
  double *Supercool = new double[n_z]; //supercooled water in soil
  double *MIce_L = new double[n_z]; //soil ice mass [mm]
  double *MLiq_L = new double[n_z]; //snow/soil liquid mass [mm]
  double *MHeat_L = new double[n_z]; //energy residual [w/m2] HM = MHeat_L
  double *MPC_L = new double[n_z]; //melting or freezing water [kg/m2] XM_L = mass of phase change

  double *SMCT_c = new double[n_z];
  double *MLiq_c = new double[n_z]; 
  double *MIce_c = new double[n_z];

  int *IndexMelt = new int[n_z]; // tracking melting/freezing index of layers
  
  int nsnow = -1;
  //compute mass of liquid/ice in soil layers in mm
  for (int i=0; i<n_z;i++) {
    if (i < nsnow) { //snow layer
      MIce_L[i] = 0; //SNICE[i];
      MLiq_L[i] = 0; //SNLIQ[i];
    }
    else {
      // MICE and MLIQ are in units of [kg/m2]
      MIce_L[i] = (SMCT[i] - SMCLiq[i]) * Dz[i] * prop.wdensity_; // [kg/m2]
      MLiq_L[i] = SMCLiq[i] * Dz[i] * prop.wdensity_;
    }
  }
  //set local variables
  
  //create copies of the current Mice and MLiq
  memcpy(MLiq_c, MLiq_L, sizeof (double) * n_z);
  memcpy(MIce_c, MIce_L, sizeof (double) * n_z);

  //Phase change between ice and liquid water
  for (int i=0; i<n_z;i++) {
    IndexMelt[i] = 0;
    SMCT_c[i] = MIce_L[i] + MLiq_L[i];
  }

  /*------------------------------------------------------------------- */
  //Soil water potential
  // SUPERCOOL is the maximum liquid water that can exist below (T - TFRZ) freezing point
  double b = -1./prop.bexp_;
  for (int i=0; i<n_z;i++) {
    if (ST[i] < prop.tfrez_) {
      double SMP = prop.lhf_ /(prop.grav_*ST[i]) * (prop.tfrez_ - ST[i]);     // [m] Soil Matrix potential
      Supercool[i] = prop.smcmax_* pow((SMP/prop.psisat_),b); //SMCMAX = porsity
      Supercool[i] = Supercool[i]*Dz[i]* prop.wdensity_; //[kg/m2];
    }
  }


  /*------------------------------------------------------------------- */
  // ****** get layer freezing/melting index ************
  for (int i=0; i<n_z;i++) {
    if (MIce_L[i] > 0 && ST[i] > prop.tfrez_) //Melting condition
      IndexMelt[i] = 1;
    else if (MLiq_L[i] > Supercool[i] && ST[i] <= prop.tfrez_)// freezing condition
      IndexMelt[i] = 2;
  }

  SoilHeatCapacity();
  /*------------------------------------------------------------------- */
  // ****** get excess or deficit of energy during phase change (use Hm) ********
  //  HC = volumetic heat capacity [J/m3/K]
  // Heat Mass = (T- Tref) * HC * DZ /Dt = K * J/(m3 * K) * m * 1/s = (J/s)*m/m3 = W/m2
  //if HeatMass < 0 --> freezing energy otherwise melting energy
  
  for (int i=0; i<n_z;i++) {
    if (IndexMelt[i] > 0) {
      MHeat_L[i] = (ST[i] - prop.tfrez_) * (HC[i] * Dz[i]) / dt;
      ST[i] = prop.tfrez_; // Note the temperature does not go below 0 until there is mixture of water and ice
    }

    if (IndexMelt[i] == 1 && MHeat_L[i] <0) {
      MHeat_L[i] = 0;
      IndexMelt[i] = 0;
    }
    
    if (IndexMelt[i] == 2 && MHeat_L[i] > 0) {
      MHeat_L[i] = 0;
      IndexMelt[i] = 0;
    }

  // compute the amount of melting or freezing water [kg/m2]. That is, how much water needs to be melted or freezed for the given energy change: MPC = MassPhaseChange
  MPC_L[i] = MHeat_L[i]*dt/prop.lhf_;
  }

  
  /*------------------------------------------------------------------- */
  // The rate of melting and freezing for snow and soil
  // mass partition between ice and water and the corresponding adjustment for the next timestep
  for (int i=0; i<n_z;i++) {
    if (IndexMelt[i] >0 && std::abs(MHeat_L[i]) >0) {
      if (MPC_L[i] >0) //melting
	MIce_L[i] = std::max(0., MIce_c[i]-MPC_L[i]);
      else if (MPC_L[i] <0) { //freezing
	if (SMCT_c[i] < Supercool[i])
	  MIce_L[i] = 0;
	else {
	  MIce_L[i] = std::min(SMCT_c[i] - Supercool[i], MIce_c[i] - MPC_L[i]);
	  MIce_L[i] = std::max(MIce_L[i],0.0);
	}
      }
    
      // compute heat residual
      // total energy available - energy consumed by phase change (ice_old - ice_new)
      double HEATR = MHeat_L[i] - prop.lhf_*(MIce_c[i]-MIce_L[i])/dt; // [W/m2] Energy Residual, last part is the energy due to change in ice mass
      MLiq_L[i] = std::max(0.,SMCT_c[i] - MIce_L[i]);

      // Temperature correction

      if (std::abs(HEATR)>0) {
	  double f = dt/(HC[i] * Dz[i]); // [m2 K/W]
	  ST[i] = ST[i] + f*HEATR; // [K] , this is computed from HeatMass = (T_n+1-T_n) * Heat_capacity * DZ/ DT
	}

	       
    }
  }
  /*
    for i in range(domain.NSNOW): #snow
    SNLIQ[i] = MLIQ_L[i] #these are already in mm, so no conversion
    SNICE[i] = MICE_L[i]
  */
  for (int i=0; i<n_z;i++) { //soil
    SMCLiq[i] =  MLiq_L[i] / (prop.wdensity_ * Dz[i]); // [-]
    SMCT[i]  = (MLiq_L[i] + MIce_L[i]) / (prop.wdensity_ * Dz[i]); // [-]
    SMCIce[i] = std::max(SMCT[i] - SMCLiq[i],0.);
  }
}

freezethaw::Properties::
Properties() :
  hcwater_ (4.188E06),
  hcice_   (2.094E06), 
  hcair_   (1004.64),
  hcsoil_  (2.00E+6),
  tcice_   (2.2), 
  lhf_     (0.3336E06),
  psisat_  (0.759),
  grav_    (9.86),
  tcwater_  (0.57),
  tcquartz_ (7.7), 
  quartz_   (0.6), 
  tcmineral_ (2.0),
  tfrez_     (273.15),
  smcmax_  (0.4), 
  bexp_  (2.9), 
  wdensity_ (1000)
{}

freezethaw::FreezeThaw::
~FreezeThaw()
{
  this->forcing_file=" ";
  //this->time = 0.;
}

#endif
