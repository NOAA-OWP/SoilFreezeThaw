#ifndef FSC_INCLUDED
#define FSC_INCLUDED

#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
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
  this->tbot = 275.15;
  this->opt_botb = 2;
  this->opt_topb = 2;
  this->forcing_file= " ";
  this->nsteps=0;
  this->ice_fraction_scheme= " ";
  this->ice_fraction_scheme_bmi = new int[1];
  this->Zd =0.0;
  this->ice_fraction_schaake =0.0;
  this->ice_fraction_xinan =0.0;
}

freezethaw::FreezeThaw::
FreezeThaw(std::string config_file)
{

  this->config_file = config_file;
  this->forcing_file= " ";
  this->lhf = 0.3336E06;
  this->ttop = 260.;
  this->tbot = 275.15;
  this->opt_botb = 2; // 1: zero thermal flux, 2: constant Temp
  this->opt_topb = 2; // 1: constant temp, 2: from a file
  
  this->ice_fraction_scheme_bmi = new int[1];
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
  this->ice_fraction_schaake =0.0;
  this->ice_fraction_xinan =0.0;
  //  if (this->ice_fraction_scheme != "BMI")
  //  GetIceFraction(); 
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

void freezethaw::FreezeThaw::
InitFromConfigFile()
{ 
  std::ifstream fp;
  fp.open(config_file);
  int n_st, n_mct, n_mcl;

  bool is_forcing_file_set = false;
  this->is_SMC_BMI_set = false;
  bool is_endtime_set = false;
  bool is_dt_set = false;
  bool is_Z_set = false;
  bool is_smcmax_set = false;
  bool is_bexp_set = false;
  bool is_quartz_set = false;
  bool is_satpsi_set = false;
  bool is_ST_set = false;
  bool is_SMCT_set = false; //total moisture content
  bool is_SMCL_set = false; //liquid moisture content
  bool is_IFS_set = false; //ice fraction scheme
    
  while (fp) {

    std::string line;
    std::string param_key, param_value, param_unit;

    std::getline(fp, line);

    int loc_eq = line.find("=") + 1;
    int loc_u = line.find("[");
    param_key = line.substr(0,line.find("="));

    bool is_unit = line.find("[") != std::string::npos;

    if (is_unit)
      param_unit = line.substr(loc_u,line.find("]")+1);
    else
      param_unit = "";

    param_value = line.substr(loc_eq,loc_u - loc_eq);
    
    if (param_key == "forcing_file") {
      this->ReadForcingData(param_value);
      is_forcing_file_set = true;
      continue;
    }
    if (param_key == "SMC_BMI") {
      this->is_SMC_BMI_set = true;
      continue;
    }
    if (param_key == "end_time") {
      this->endtime = std::stod(param_value);

      if (param_unit == "[d]" || param_unit == "[day]") 
	this->endtime *= 86400;
      else if (param_unit == "[s]" || param_unit == "[sec]")
	this->endtime *= 1.0;
      else if (param_unit == "[h]" || param_unit == "[hr]" || param_unit == "") // defalut time unit is hour
	this->endtime *= 3600.0;

      is_endtime_set = true;
      continue;
    }
    if (param_key == "dt") {
      this->dt = std::stod(param_value);
      if (param_unit == "[d]" || param_unit == "[day]")
	this->dt *= 86400;
      else if (param_unit == "[s]" || param_unit == "[sec]")
	this->dt *= 1.0;
      else if (param_unit == "[h]" || param_unit == "[hr]" || param_unit == "") // defalut time unit is hour
	this->dt *= 3600.0;
      
      is_dt_set = true;
      continue;
    }
    if (param_key == "Z") {
      std::vector<double> vec = ReadVectorData(param_value);
      
      this->Z = new double[vec.size()];
      for (unsigned int i=0; i < vec.size(); i++)
	this->Z[i] = vec[i];
      this->nz = vec.size();
      this->Zd = this->Z[this->nz-1];
      is_Z_set = true;
      continue;
    }
    if (param_key == "soil_params.smcmax") {
      this->smcmax = std::stod(param_value);
      is_smcmax_set = true;
      continue;
    }
    if (param_key == "soil_params.b") {
      this->bexp = std::stod(param_value);
      std::string bexp_unit = line.substr(loc_u+1,line.length());
      assert (this->bexp > 0);
      is_bexp_set = true;
      continue;
    }
    if (param_key == "soil_params.quartz") {
      this->quartz = std::stod(param_value);
      assert (this->quartz > 0);
      is_quartz_set = true;
      continue;
    }
    if (param_key == "soil_params.satpsi") {  //Soil saturated matrix potential
      this->satpsi = std::stod(param_value);
      is_satpsi_set = true;
      continue;
    }
    if (param_key == "soil_temperature") {
      std::vector<double> vec = ReadVectorData(param_value);
      this->ST = new double[vec.size()];
      for (unsigned int i=0; i < vec.size(); i++)
	this->ST[i] = vec[i];
      n_st = vec.size();
      
      is_ST_set = true;
      continue;

    }
    if (param_key == "soil_total_moisture_content") {
      std::vector<double> vec = ReadVectorData(param_value);
      this->SMCT = new double[vec.size()];
      for (unsigned int i=0; i < vec.size(); i++)
	this->SMCT[i] = vec[i];
      n_mct = vec.size();
      is_SMCT_set = true;
      continue;
    }
    if (param_key == "soil_liquid_moisture_content") {
      std::vector<double> vec = ReadVectorData(param_value);
      this->SMCLiq = new double[vec.size()];
      for (unsigned int i=0; i < vec.size(); i++) {
	//	assert (this->SMCT[i] >= vec[i]);
	this->SMCLiq[i] = vec[i];
      }
      n_mcl = vec.size();
      is_SMCL_set = true;
      continue;
    }
    if (param_key == "ice_fraction_scheme") {
      this->ice_fraction_scheme = param_value;
      is_IFS_set = true;
      continue;
    }
  }
  
  fp.close();
  
  // simply allocate space for SMCLiq and SMCT arrays, as they will be set through CFE_BMI
  if (this->is_SMC_BMI_set && is_Z_set) {
    this->SMCT = new double[this->nz]();
    this->SMCLiq = new double[this->nz]();
    n_mct = this->nz;
    n_mcl = this->nz;
    is_SMCT_set = true;
  }
  
  if (!is_forcing_file_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("Forcing file not set in the config file!");
  }
  if (!is_endtime_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("End time not set in the config file!");
  }

  if (!is_dt_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("Time step (dt) not set in the config file!");
  }
  if (!is_Z_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("Z not set in the config file!");
  }
  if (!is_smcmax_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("smcmax not set in the config file!");
  }
  if (!is_bexp_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("bexp (Clapp-Hornberger's parameter) not set in the config file!");
  }
  if (!is_quartz_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("quartz (soil parameter) not set in the config file!");
  }
  if (!is_satpsi_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("satpsi not set in the config file!");
  }
  if (!is_ST_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("Soil temperature not set in the config file!");
  }
  if (!is_SMCT_set && !this->is_SMC_BMI_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("Total soil moisture content not set in the config file!");
  }
  if (!is_SMCL_set && !this->is_SMC_BMI_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("Liquid soil moisture content not set in the config file!");
  }
  if (!is_IFS_set) {
    std::cout<<"Config file: "<<this->config_file<<"\n";
    throw std::runtime_error("Ice fraction scheme not set in the config file!");
  }

  // check if the size of the input data is consistent
  assert (n_st == this->nz);
  assert (n_mct == this->nz);
  assert (n_mcl == this->nz);
}


std::vector<double> freezethaw::FreezeThaw::
ReadVectorData(std::string key)
{
  int pos =0;
  std::string delimiter = ",";
  std::vector<double> value(0);
  std::string z1 = key;

  while (z1.find(delimiter) != std::string::npos) {
    pos = z1.find(delimiter);
    //std::ostringstream z_vv;
    //z_vv << std::setprecision(8) << z1.substr(0, pos);
    std::string z_v = z1.substr(0, pos);

    value.push_back(stod(z_v.c_str()));

    z1.erase(0, pos + delimiter.length());
    if (z1.find(delimiter) == std::string::npos)
      value.push_back(stod(z1));
  }

  return value;
}

void freezethaw::FreezeThaw::
ReadForcingData(std::string forcing_file)
{
  std::ifstream fp;
  fp.open(forcing_file);
  if (!fp) {
    cout<<"file "<<forcing_file<<" doesn't exist. \n";
    abort();
  }
  
  std::vector<double> Time_v(0.0);
  std::vector<double> GT_v(0.0);
  std::vector<string> vars;
  std::string line, cell;
  
  //read first line of strings which contains forcing variables names.
  std::getline(fp, line);
  std::stringstream lineStream(line);
  int ground_temp_index=-1;
  
  while(std::getline(lineStream,cell, ',')) {
    vars.push_back(cell);
  }

  for (unsigned int i=0; i<vars.size();i++) {
    if (vars[i] ==  "TMP_ground_surface")
      ground_temp_index = i;
  }

  if (ground_temp_index <0)
    ground_temp_index = 6; // 6 is the air temperature column, if not coupled and ground temperatgure is not provided
    
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

      if (count % len_v == ground_temp_index) {
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

  // this is needed to make sure external calls (such as CFE BMI) don't exceed the length of the SFT forcing data
  this->total_nsteps = size_v;
}

void freezethaw::FreezeThaw::
GetIceFraction()
{

  double val = 0;
  this->ice_fraction_schaake = 0.0; // set it to zero
  this->ice_fraction_xinan = 0.0;
  
  if (this->ice_fraction_scheme == "Schaake") {
    *this->ice_fraction_scheme_bmi = 1;
  }
  else if (this->ice_fraction_scheme == "Xinanjiang") {
    *this->ice_fraction_scheme_bmi = 2;
  }
  
  if (*this->ice_fraction_scheme_bmi == SurfaceRunoffScheme::Schaake) {
    val = this->SMCIce[0]*this->Z[0];
    for (int i =1; i < nz; i++) {
      val += this->SMCIce[i] * (this->Z[i] - this->Z[i-1]);
    }
    assert (this->ice_fraction_schaake <= this->Zd);
    this->ice_fraction_schaake = val;
  }
 else if (*this->ice_fraction_scheme_bmi == SurfaceRunoffScheme::Xinanjiang) {
    double fice = std::min(1.0, this->SMCIce[0]/this->smcmax);
    double A = 4.0; // taken from NWM SOILWATER subroutine
    double fcr = std::max(0.0, std::exp(-A*(1.0-fice)) - std::exp(-A)) / (1.0 - std::exp(-A));
    this->ice_fraction_xinan = fcr;
  }
  else {
    throw std::runtime_error("Ice Fraction Scheme not specified either in the config file nor set by CFE BMI. Options: Schaake or Xinanjiang!");
    }
}
  
double freezethaw::FreezeThaw::
GetDt()
{
  return this->dt;
}

void freezethaw::FreezeThaw::
Advance()
{
  assert (this->nsteps < this->total_nsteps);
  
  if (this->is_SMC_BMI_set) {
    for (int i=0; i<this->nz;i++) {
      this->SMCLiq[i] = this->SMCT[i];
      this->SMCIce[i] = 0.0;
    }
    this->is_SMC_BMI_set = false;
  }

  /*
  if (this->is_SMC_BMI_set) {
    for (int i=0; i<this->nz;i++) {
      this->SMCLiq[i] = this->SMCT[i] - this->SMCIce[i];
      //      this->SMCIce[i] = 0.0;
    }
  }*/
  
  // Update Thermal conductivities, note the soil heat flux update happens in the PhaseChange module, so no need to update here
  ThermalConductivity(); // initialize thermal conductivities

  // Update volumetric heat capacity
  SoilHeatCapacity();
  
  // call Diffusion Eq. first to get the updated soil temperatures
  SolveDiffusionEq();

  // call phase change to get water/ice partition for the updated temperature
  PhaseChange();

  this->time += this->dt;

  GetIceFraction();

  
  assert (this->ST[0] >150.0); // getting temperature below 200 would mean the space resolution is too fine and time resolution is too coarse
  
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
    assert (this->Z[0] >0);
    double ghf = - TC[0] * (surfT  - GT[this->nsteps]) / this->Z[0];  //temperature from a file
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
    const int nz = this->shape[0];

    // local variables
    std::vector<double> Flux(nz);
    std::vector<double> AI(nz);
    std::vector<double> BI(nz);
    std::vector<double> CI(nz);
    std::vector<double> RHS(nz);
    std::vector<double> Lambd(nz);
    std::vector<double> X(nz);
    double botflux=0.0;
    double h1 =0.0, h2 =0.0;
    
    // compute matrix coefficient using Crank-Nicolson discretization scheme
    for (int i=0;i<nz; i++) {
      if (i == 0) {
	h1 = Z[i];
	double dtdz  = (ST[i+1] - ST[i])/ h1;
	Lambd[i] = dt/(2.0 * h1 * HC[i]);
	double ghf = this->GroundHeatFlux(ST[i]);
	Flux[i] = Lambd[i] * (TC[i] * dtdz + ghf);
      }
      else if (i < nz-1) {
	h1 = Z[i] - Z[i-1];
        h2 = Z[i+1] - Z[i];
	Lambd[i] = dt/(2.0 * h2 * HC[i]);
	double a_ = - Lambd[i] * TC[i-1] / h1;
        double c_ = - Lambd[i] * TC[i] / h2;
	double b_ = 1 + a_ + c_;
	Flux[i] = -a_ * ST[i-1] + b_ * ST[i] - c_ * ST[i+1];
      }
      else if (i == nz-1) {
	h1 = Z[i] - Z[i-1];
	Lambd[i] = dt/(2.0 * h1 * HC[i]);
	if (this->opt_botb == 1) 
	  botflux = 0.;
	else if (this->opt_botb == 2) {
	  double dtdz1 = (ST[i] - tbot) / h1;
	  botflux  = - TC[i] * dtdz1;
	}
	double dtdz = (ST[i] - ST[i-1] )/ h1;
	Flux[i]  = Lambd[i] * (-TC[i]*dtdz  + botflux);
      }
    }
    
    // put coefficients in the corresponding vectors A,B,C, RHS
    for (int i=0; i<nz;i++) {
      if (i == 0) {
	AI[i] = 0;
	CI[i] = - Lambd[i] *TC[i]/Z[i];
	BI[i] = 1 - CI[i];
      }
      else if (i < nz-1) {
	AI[i] = - Lambd[i] * TC[i-1]/(Z[i] - Z[i-1]);
	CI[i] = - Lambd[i] * TC[i]/(Z[i+1] - Z[i]);
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
    for (int i=0; i<nz;i++) {
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

  double tcmineral = this->quartz > 0.2 ? 2.0 : 3.0; //TC of other mineral
  double tcquartz = 7.7; // TC of Quartz
  double tcwater  = 0.57; // TC of water
  double tcice    = 2.2;  // thermal conductiviyt of ice
  
  for (int i=0; i<n_z;i++) {
    double sat_ratio = SMCT[i]/ this->smcmax;

    //TC of solids Eq. (10) Peters-Lidard
    double tc_solid = pow(tcquartz,this->quartz) * pow(tcmineral, (1. - this->quartz));

    //SATURATED THERMAL CONDUCTIVITY
    
    //UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
    double x_unfrozen= 1.0; //prevents zero division
    if (this->SMCT[i] > 0)
      x_unfrozen = this->SMCLiq[i] / this->SMCT[i]; // (phi * Sliq) / (phi * sliq + phi * sice) = sliq/(sliq+sice) 
    
    double xu = x_unfrozen * this->smcmax; // unfrozen volume fraction
    double tc_sat = pow(tc_solid,(1. - this->smcmax)) * pow(tcice, (this->smcmax - xu)) * pow(tcwater,xu);
    
    //DRY THERMAL CONDUCTIVITY
    double gammd = (1. - this->smcmax)*2700.; // dry density
    double tc_dry = (0.135* gammd+ 64.7)/ (2700. - 0.947* gammd);
    
    // Kersten Number
    
    double KN;
    if ( (SMCLiq[i] + 0.0005) < SMCT[i])
      KN = sat_ratio; // for frozen soil
    else {
      if (sat_ratio > 0.1)
	KN = log10(sat_ratio) + 1.;
      else if (sat_ratio > 0.05)
	KN = 0.7 * log10(sat_ratio) + 1.;
      else
	KN = 0.0;
    }
    
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
    HC[i] = SMCLiq[i]*prop.hcwater_ + sice*prop.hcice_ + (1.0-this->smcmax)*prop.hcsoil_ + (this->smcmax-SMCT[i])*prop.hcair_;
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
  double lam = -1./(this->bexp);
  for (int i=0; i<n_z;i++) {
    if (ST[i] < prop.tfrez_) {
      double SMP = prop.lhf_ /(prop.grav_*ST[i]) * (prop.tfrez_ - ST[i]);     // [m] Soil Matrix potential
      Supercool[i] = this->smcmax* pow((SMP/this->satpsi), lam); //SMCMAX = porsity
      Supercool[i] = Supercool[i]*Dz[i]* prop.wdensity_; //[kg/m2];
    }
  }


  /*------------------------------------------------------------------- */
  // ****** get layer freezing/melting index ************
  for (int i=0; i<n_z;i++) {
    if (MIce_L[i] > 0 && ST[i] > prop.tfrez_) //Melting condition
      IndexMelt[i] = 1;
    else if (MLiq_L[i] > Supercool[i] && ST[i] <= prop.tfrez_)// freezing condition in NoahMP
      IndexMelt[i] = 2;
  }

  //SoilHeatCapacity();
  
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
  //  tcice_   (2.2), 
  lhf_     (0.3336E06),
  //  satpsi_  (0.759),
  grav_    (9.86),
  //tcwater_  (0.57),
  //tcquartz_ (7.7), 
  //  quartz_   (0.6), 
  //  tcmineral_ (2.0),
  tfrez_     (273.15),
  //  bexp_  (2.9), 
  wdensity_ (1000)
{}

freezethaw::FreezeThaw::
~FreezeThaw()
{
  this->forcing_file=" ";
  //this->time = 0.;
}

#endif
