#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <iostream>
#include <fstream>

#include "../bmi/bmi.hxx"
#include "../include/bmi_soil_freeze_thaw.hxx"
#include "../include/soil_freeze_thaw.hxx"
#include <cmath>



std::vector<double> ReadForcingData(std::string config_file);

/************************************************************************
   The code simulates a standalone run of the Soil Freeze-thaw model.
   Benchmark: Comparison of the ice fraction is made with the existing (already ran) golden test
************************************************************************/


int main(int argc, const char *argv[])
{

  /************************************************************************
      A configuration file is required for running this model through BMI
  ************************************************************************/
  if(argc<=1){
    printf("make sure to include a path to the SFT config file\n");
    exit(1);
  }

  /************************************************************************
  Creating SFT object
  ************************************************************************/
  BmiSoilFreezeThaw ftm_bmi_model;
 
  /************************************************************************
      Initializing the BMI model for Soil freeze-thaw model
  ************************************************************************/
  
  const char *cfg_file_ftm = argv[1];
  ftm_bmi_model.Initialize(cfg_file_ftm);

  //Read ground temperature data for SFT
  std::vector<double> ground_temp = ReadForcingData(cfg_file_ftm);
  
  /************************************************************************
    Now loop through time and call the models with the intermediate get/set
  ************************************************************************/

  // get time steps
  double endtime = ftm_bmi_model.GetEndTime();
  double timestep = ftm_bmi_model.GetTimeStep ();
  int nsteps = int(endtime/timestep); // total number of time steps
  //  nsteps = 24566;
  
  double *ice_fraction = new double[nsteps];
  std::vector<double> ice_fraction_golden;
  double ice_frac_v = 0;

  bool golden_test = false;
  std::ofstream outfile;
 
  std::string filename = "./test/file_golden.csv";
  
  if (golden_test) {
    
    outfile.open(filename, fstream::out);
    
    if (!outfile)
      std::cout<<"error"<<"\n";

    outfile << "Time [h],ice_fraction" << "\n"; 
  }
  else {
    
    std::ifstream infile;
    infile.open(filename);
    
    if (!infile)
      std::cout<<"Can't open the file "<< filename<<"\n";

    int a;
    double b;
    char c;
    std::string line;
    //get the header
    std::getline(infile,line);

    while ((infile >> a >> c >> b) && (c == ',')) {
      ice_fraction_golden.push_back(b);
    }
    
    
    infile.close();
  }
  
  for (int i = 0; i < nsteps; i++) {
    
    ftm_bmi_model.SetValue("ground_temperature", &ground_temp[i]);

    ftm_bmi_model.Update(); // Update model

    ftm_bmi_model.GetValue("ice_fraction_schaake",&ice_frac_v);

    ice_fraction[i] = ice_frac_v;
    
    if (golden_test)
      outfile << i+1 << "," <<ice_frac_v << "\n";
    
  }

  outfile.close();
  
  /*********************** Comnpare against golden test ********************************/

  if (!golden_test) {
    bool test_status = true;
    
    double  err_frozen_frac_mm = 0;
    for (int i=0; i<nsteps;i++) {
      err_frozen_frac_mm += round(fabs(ice_fraction_golden[i] - ice_fraction[i]) *100000.)/100000.; //truncate the error at 5 decimal places
    }
    
    //RMSE
    err_frozen_frac_mm  = std::pow(err_frozen_frac_mm/nsteps,0.5);
    
    if (err_frozen_frac_mm < 1.e-2)
      test_status &= true;
    else
      test_status &= false;
    
    std::string passed = test_status > 0 ? "Yes" : "No";
    
    
    std::cout<<"*********************************************************\n";
    std::cout<<" Test passed = "<<passed<<" \n Frozen fraction error = "<<err_frozen_frac_mm<<"\n";
    std::cout<<"*********************************************************\n";
  }
  else {
    std::cout<<"Golden test created... see "<<filename<<"\n";
  }
  
  /************************************************************************
    Finalize SFT BMI model
  ************************************************************************/
  
  ftm_bmi_model.Finalize();
  
  return 0;
}


std::vector<double>
ReadForcingData(std::string config_file)
{
  // get the forcing file from the config file

  std::ifstream file;
  file.open(config_file);

  if (!file) {
    std::stringstream errMsg;
    errMsg << config_file << " does not exist";
    throw std::runtime_error(errMsg.str());
  }

  std::string forcing_file;
  bool is_forcing_file_set=false;
  
  while (file) {
    std::string line;
    std::string param_key, param_value;

    std::getline(file, line);

    int loc_eq = line.find("=") + 1;
    param_key = line.substr(0, line.find("="));
    param_value = line.substr(loc_eq,line.length());

    if (param_key == "forcing_file") {
      forcing_file = param_value;
      is_forcing_file_set = true;
      break;
    }
  }

  if (!is_forcing_file_set) {
    std::stringstream errMsg;
    errMsg << config_file << " does not provide forcing_file";
    throw std::runtime_error(errMsg.str());
  }
  
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

  return GT_v;
 
}
