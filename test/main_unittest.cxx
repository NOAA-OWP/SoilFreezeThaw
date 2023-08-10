/*
  @author Ahmad Jan (ahmad.jan@noaa.gov)
  Includes unit test for bmi components and run model for 2 days with hourly timestep to test the computed ice fraction using schaake runoff scheme
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <iomanip>      // std::setprecision
#include "../bmi/bmi.hxx"
#include "../include/bmi_soil_freeze_thaw.hxx"
#include "../include/soil_freeze_thaw.hxx"

#define FAILURE 0
#define VERBOSITY 1

#define GREEN "\033[32m"
#define RED   "\033[31m"
#define BLUE  "\033[34m"
#define RESET "\033[0m"


int main(int argc, char *argv[])
{
  BmiSoilFreezeThaw model, model_cyc, model_calib;

  if (argc != 2) {
    printf("Usage: ./run_unittest.sh \n\n");
    printf("Run the frozensoilcxx model through its BMI with a configuration file.\n");
    return FAILURE;
  }

  std::cout<<"\n**************** BEGIN SoilFreezeThaw BMI UNIT TEST *******************\n";

  model.Initialize(argv[1]);
  model_cyc.Initialize(argv[1]);
  model_calib.Initialize(argv[1]);

  std::cout<<"\n**************** TEST VALUES ************************************\n";
  
  int    nz          = 4;
  double endtime     = 86400; //1035.91*86400.;
  double timestep    = 3600;
  bool   test_status = true;
  
  std::vector<string> bmi_input_vars = {"ground_temperature", "soil_moisture_profile"};
  std::vector<string> bmi_output_vars = {"ice_fraction_schaake", "ice_fraction_xinan", "num_cells",
					 "soil_temperature_profile", "soil_ice_fraction", "ground_heat_flux"};

  int num_input_vars  = bmi_input_vars.size();
  int num_output_vars = bmi_output_vars.size();
  
  int nbytes_input[] = {sizeof(double), int(nz * sizeof(double))};
  int nbytes_output[] = {sizeof(double), sizeof(double), sizeof(int), int(sizeof(double) * nz), sizeof(double),
			 sizeof(double)};
  
  double soil_moisture_profile[] = {0.389,0.396,0.397,0.397}; // total_moisture_content
  double soil_T[] = {280.15,280.15,280.15,280.15}; //soil temperature
  
  std::cout<<"Num cells:           "<<nz<<"\n";
  std::cout<<"End time:            "<<endtime<<"\n";
  std::cout<<"Num input vars:      "<<num_input_vars<<"\n";
  std::cout<<"Num output vars:     "<<num_output_vars<<"\n";

  std::cout<<"\nPulling information from BMI\n************************************\n";
  std::string model_name;
  int count_in = 0;
  int count_out = 0;
  std::vector<std::string> names_in;
  std::vector<std::string> names_out;
  

  // Test get_component_name()
  model_name = model.GetComponentName();
  if (VERBOSITY)
    std::cout<<"Model name: "<< model_name <<"\n";

  // Test GetInputItemCount
  count_in = model.GetInputItemCount();
  if (VERBOSITY)
    std::cout<<"Input item count: "<< count_in<<"\n";
  if (count_in == num_input_vars)
    test_status &= true;
  else {
    test_status &= false;
    std::string passed = test_status == true ? "Yes" : "No";
    std::cout<<"Test passed: "<<passed<<"\n";
    std::stringstream errMsg;
    errMsg << "Number of input variables are different. "<< count_in << " != "<< num_input_vars << "\n";
    throw std::runtime_error(errMsg.str());
  }

  // Test GetInputVarNames 
  names_in = model.GetInputVarNames();
  if (VERBOSITY) {
    std::cout<<"Input variable names \n";
    for (int i=0; i<count_in; i++)
      std::cout<<i<<" "<<names_in[i]<<"\n";
  }

  std::cout<<"**************************************** \n";
  // Test GetOutputItemCount
  count_out = model.GetOutputItemCount();
  if (VERBOSITY)
    std::cout<<"Output item count: "<< count_out<<"\n";

  // Test GetOutputVarNames
  names_out = model.GetOutputVarNames();
  if (VERBOSITY) {
    std::cout<<"Output variable names "<<names_out.size()<<"\n";
    for (int i=0; i<count_out; i++)
      std::cout<<i<<" "<<names_out[i]<<"\n";
  }
  if (count_out == num_output_vars)
    test_status &= true;
  else {
    test_status &= false;
    std::string passed = test_status == true ? "Yes" : "No";
    std::cout<<"Test passed: "<<passed<<"\n";
    std::stringstream errMsg;
    errMsg << "Number of output variables are different. "<< count_out <<" != "<< num_output_vars <<"\n";
    throw std::runtime_error(errMsg.str());
  }
  
  // Test BMI: VARIABLE INFORMATION FUNCTIONS
  std::cout<<"\n**************** TEST BMI VARIABLE INFORMATION FUNCTIONS\n***************************\n";
  int grid, itemsize, nbytes;
  std::string location;
  std::string units;
  std::string vartype;
  
  // Loop through both input and output variables and call GetVar* functions
  for (int i=0; i<count_in; i++) {
    std::string var_name = names_in[i];
    if (VERBOSITY)
      std::cout<<"Input var_name: "<< var_name <<"\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_grid()
    grid = model.GetVarGrid(var_name);
    if (VERBOSITY)
      std::cout<<"Grid: "<< grid <<"\n";

    if (grid >=0)
      test_status &= true;
    else {
      test_status &= false;
      std::string passed = test_status == true ? "Yes" : "No";
      std::cout<<"Test passed: "<<passed<<"\n";
      std::stringstream errMsg;
      errMsg << "grid < 0 \n";
      throw std::runtime_error(errMsg.str());
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_itemsize()
    itemsize = model.GetVarItemsize(var_name);
    if (VERBOSITY)
      std::cout<<"Itemsize: "<< itemsize <<"\n";
    if (itemsize >0)
      test_status &= true;
    else {
      test_status &= false;
      std::string passed = test_status == true ? "Yes" : "No";
      std::cout<<"Test passed: "<<passed<<"\n";
      std::stringstream errMsg;
      errMsg << "itemsize < 0 \n";
      throw std::runtime_error(errMsg.str());
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_location()
    location = model.GetVarLocation(var_name);
    if ( location == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" location: "<< location<<"\n";
    if (location == "")
      test_status &= false;
    else
      test_status &= true;
	
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_units()
    units = model.GetVarUnits(var_name);
    if (VERBOSITY)
      std::cout<<" units: ["<< units <<"]\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_type()
    vartype = model.GetVarType(var_name);
    if (VERBOSITY)
      std::cout<<" type: "<< vartype <<"\n";
    if (location == "")
      test_status &= false;
    else
      test_status &= true;
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // get_var_nbytes()
    nbytes = model.GetVarNbytes(var_name);
    if (nbytes == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<" nbytes: "<< nbytes <<"\n";

    if (var_name == bmi_input_vars[i]) {
      if (nbytes == nbytes_input[i])
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::stringstream errMsg;
	errMsg << "Number of bytes for input var"<<var_name<< " should be "<<nbytes_input[i]<<"\n";
	throw std::runtime_error(errMsg.str());
      }
    }
    else {
      std::stringstream errMsg;
      errMsg << "Variable name"<< var_name<<" should be: soil_moisture_profile or ground_temperature \n";
      throw std::runtime_error(errMsg.str());

    }
  }

  if (VERBOSITY)
    std::cout<<"\n*****************************************\n";
  
  for (int i=0; i<count_out; i++) {
    std::string var_name = names_out[i];
    if (VERBOSITY)
      std::cout<<"Output var_name: "<< var_name <<"\n";
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_grid() 
    grid = model.GetVarGrid(var_name);
    if (grid == -1) return -1;
    if (VERBOSITY)
      std::cout<<"Grid: "<< grid <<"\n";

    if (grid >=0)
      test_status &= true;
    else
      test_status &= false;
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_itemsize()
    itemsize = model.GetVarItemsize(var_name);
    if (itemsize == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<"Itemsize: "<< itemsize <<"\n";

    if (itemsize >0)
      test_status &= true;
    else
      test_status &= false;
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_location()
    location = model.GetVarLocation(var_name);
    if ( location == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" location:"<< location<<"\n";

    if (location == "")
      test_status &= false;
    else
      test_status &= true;
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_units()
    units = model.GetVarUnits(var_name);
    if (VERBOSITY)
      std::cout<<" units: ["<< units <<"]\n";
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_type()
    vartype = model.GetVarType(var_name);
    if (vartype == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" type: "<< vartype <<"\n";
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // get_var_nbytes()
    nbytes = model.GetVarNbytes(var_name);
    if (nbytes == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<" nbytes: "<< nbytes<<"\n";

    if (var_name == bmi_output_vars [i]) {
      if (nbytes == nbytes_output[i])
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::stringstream errMsg;
	errMsg << "Number of bytes for output var"<<var_name<< " should be "<<nbytes_output[i]<<"\n";
	throw std::runtime_error(errMsg.str());
      }
    }
    else {
      std::stringstream errMsg;
      errMsg << "Variable name"<< var_name<<" should be: ice_fraction_schaake or ice_fraction_xinan or num_cell \n";
      throw std::runtime_error(errMsg.str());

    }
    
  }

  // Test BMI: MODEL GRID FUNCTIONS
  std::cout<<"\n \n**************** TEST BMI GRID FUNCTIONS***********************\n";
  int grid_id[] = {0,1,2};
  int grid_size_test[] = {1,1,nz};
  int grid_rank, grid_size;
  std::string grid_type;

  for (int i=0; i< 3; i++) {
    if (VERBOSITY)  
      std::cout<<"Grid id "<< grid_id[i] <<"\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_grid_rank()
    grid_rank = model.GetGridRank(grid_id[i]);
    if (grid_rank == FAILURE) return FAILURE;
    if (VERBOSITY)
      std::cout<<" rank: "<<grid_rank<<"\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_grid_size
    grid_size = model.GetGridSize(grid_id[i]);
    if (grid_size == grid_size_test[i]) {
      test_status &= true;
      if (VERBOSITY)
	std::cout<<" grid size: "<<grid_size<<"\n";
    }
    else {
      test_status &= false;
      std::string passed = test_status == true ? "Yes" : "No";
      std::cout<<"Test passed: "<<passed<<"\n";
      std::stringstream errMsg;
      errMsg << "Grid size of should be "<<nz<<"\n";
      throw std::runtime_error(errMsg.str());
    }
  }

  // Test BMI: TIME FUNCTIONS
  std::cout<<"\n \n**************** TEST BMI TIME FUNCTIONS***********************\n";
  double time = 0.0;
  double dt = 0.0;
  std::string units_time;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Test get_start_time()
  time = model.GetStartTime();
  if (VERBOSITY)
    std::cout<<" start time: "<< time<<"\n";

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Test get_end_time()
  time = model.GetEndTime();
  if (time == endtime) {
    test_status &= true;
    if (VERBOSITY)
      std::cout<<" end time [h]: "<< time /3600.<<"\n";
  }
  else {
    test_status &= false;
    std::string passed = test_status == true ? "Yes" : "No";
    std::cout<<"Test passed: "<<passed<<"\n";
    std::stringstream errMsg;
    errMsg << "End time should be ["<<endtime/3600.<<" hours] \n";
    throw std::runtime_error(errMsg.str());
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Test get_time_step()
  dt = model.GetTimeStep();
  if (dt == timestep) {
    if (VERBOSITY)
      std::cout<<" timestep: "<< dt<<"\n";
  }
  else
    test_status &= false;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Test get_time_units()
  units_time = model.GetTimeUnits();
  if (units_time == "s") {
    test_status &= true;
    if (VERBOSITY)
      std::cout<<" time units: "<< units_time <<"\n";
  }
  else {
    test_status &= false;
    std::string passed = test_status == true ? "Yes" : "No";
    std::cout<<"Test passed: "<<passed<<"\n";
    std::stringstream errMsg;
    errMsg << "Time units should be seconds [s], but the model returned"<<units_time<<"\n";
    throw std::runtime_error(errMsg.str());
  }

  std::cout<<GREEN<<"\n";
  std::string passed = test_status > 0 ? "Yes" : "No";
  std::cout<<"\n| *************************************** \n";
  std::cout<<"All tests passed at this point: "<<passed<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Test get_current_time()
  time = model.GetCurrentTime();
  if (VERBOSITY)
    std::cout<<" current time: "<< time<<"\n";

  // Test BMI: GET VALUE FUNCTIONS
  std::cout<<"\n\n************** TEST BMI GETTER SETTER FUNCTIONS *************************\n";
  
  std::cout<<"********** Input variables ***************** \n";
  // Loop through both input and output variables and call get/set_value_*()
  
  for (int i=0; i<count_in; i++) {
    std::string var_name = names_in[i];
    std::cout<<"Variable name: "<< var_name <<"\n";
    
    double *var = new double[nz];
    double *dest = new double[nz];
    int indices[] = {0,1,2,3};
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_value() at each timestep
    model.GetValue(var_name, &(var[0]));
    std::cout<<" Get value: "<< var[0] <<"\n";
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_value_at_indices()
    model.GetValueAtIndices(var_name, dest, indices, nz);
    std::cout<<" Get value at indices: " << dest[0]<<"\n";
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_value_ptr()
    double *var_ptr = new double[nz];
    var_ptr= (double*) model.GetValuePtr(var_name);
    std::cout<<" Get value ptr: "<<*var_ptr<<"\n";
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // compare benchmark value with bmi GetValue
    if (var_name == "soil_moisture_profile") {
      bool check = false;
      for (int i1 =0; i1 < nz; i1++) {
	assert (fabs(var[i] - soil_moisture_profile[i]) < 0.001);
	check=true;
      }
      if (check)
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::stringstream errMsg;
	errMsg << "Soil moisture should be: "<<soil_moisture_profile[0]<<" "<<soil_moisture_profile[1]<<" "
	       <<soil_moisture_profile[2]<<" "<<soil_moisture_profile[3]<<"\n";
	throw std::runtime_error(errMsg.str());
      }
    }
    else if (var_name == "ground_temperature") {
    // Go ahead and test set_value_*() for last time step here
    // Test BMI set_value()
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    double var_new = 286.6;
    double var_new_up = 0.0;
    model.SetValue(var_name, &(var_new));
    
    std::cout<<" Set value: "<< var_new <<"\n";
    // get_value to see if changed
    model.GetValue(var_name, &var_new_up);
    std::cout<<" Get value: "<< var_new_up <<"\n";    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test BMI set_value_at_indices()
    double dest_new = 281.3;
    double dest_new_up = 0.0;
    
    model.SetValueAtIndices(var_name, &indices[0], nz, &dest_new);
    
    std::cout<<" Set value at indices: "<< dest_new <<"\n";
    // get_value_at_indices to see if changed
    model.GetValueAtIndices(var_name, &dest_new_up,  &indices[0], nz);
    std::cout<<" Get value at indices: "<< dest_new_up <<"\n";
    
    if (dest_new == dest_new_up)
      test_status &= true;
    else
      test_status &= false;

    }
  }

  passed = test_status > 0 ? "Yes" : "No";
  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| All tests passed at this point: "<<passed<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";
  
  std::cout<<"************* Output variables ***************** \n";
  
  for (int i=0; i<count_out; i++) {
    
    std::string var_name = names_out[i];
    std::cout<<"variable name: "<< var_name <<" "<<test_status<<"\n";
    
    if (var_name.compare("num_cells") == 0) {
      int nz_v = 0;
      model.GetValue(var_name, &nz_v);
      std::cout<<" Get value: "<< nz_v <<"\n";
      if (nz == nz_v)
	test_status &= true;
      else
	test_status &= false;
    }
    else if (var_name.compare("ice_fraction_schaake") == 0 || var_name.compare("ice_fraction_xinan") == 0
	     || var_name.compare("soil_ice_fraction") == 0) {
      double val_v = 0;
      double val = 0.0; // benchmark value
      model.GetValue(var_name, &val_v);
      std::cout<<" Get value: "<< val_v <<"\n";
      if (val == val_v)
	test_status &= true;
      else
	test_status &= false;
    }
    else if (var_name.compare("soil_temperature_profile") == 0) {
      double *var = new double[nz];
      model.GetValue(var_name, &(var[0]));
      bool check = false;
      for (int i1 =0; i1 < nz; i1++) {
	assert (fabs(var[i] - soil_T[i]) < 0.001);
	check=true;
      }
      if (check)
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::stringstream errMsg;
	errMsg << "Soil temperature should be: "<<soil_T[0]<<" "<<soil_T[1]<<" "
	       <<soil_T[2]<<" "<<soil_T[3]<<"\n";
	throw std::runtime_error(errMsg.str());
      }
    }
    
  }

  
  passed = test_status > 0 ? "Yes" : "No";
  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| All BMI Tests passed: "<<passed<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";


  // Run the model for 480 timesteps, i.e., two days and compare soil_ice_fraction with the known ice_fraction
  
  int nstep = 480;
  std::cout<<"*********** Timestepping now.......\n";
  std::cout<<"Total timesteps = "<< nstep <<" (20 days)\n";

  double ground_temp = 280.15;
  for (int n=0; n<nstep; n++) {
    
    if (n < 200)
      ground_temp -= 0.5;
    else
      ground_temp += 0.5;
      
    std::cout<<"------------------------------------------------------ \n";
    std::cout<<"Timestep | "<< n <<", ground temp = "<< ground_temp <<"\n";
    std::cout<<"------------------------------------------------------ \n";
    model_cyc.SetValue("ground_temperature", &ground_temp);
    model_cyc.Update();    
  }
  
  double total_ice_content = 0.30072303607;  // benchmark value  in meters (this is computed based on the runoff scheme)
  double soil_ice_fraction = 0.379269814694; // benchmark value [-]
  double total_ice_content_sim, soil_ice_fraction_sim;
  
  model_cyc.GetValue("ice_fraction_schaake",&total_ice_content_sim);

  model_cyc.GetValue("soil_ice_fraction",&soil_ice_fraction_sim);
  
  double err_ice_content_mm  = fabs(total_ice_content - total_ice_content_sim) *1000;
  double err_ice_fraction = fabs(soil_ice_fraction - soil_ice_fraction_sim);
    
  if (err_ice_content_mm < 1E-3 || err_ice_fraction < 1E-3)
    test_status &= true;
  else
    test_status &= false;

  std::cout<<"\n*********************************************************\n";
  std::cout<<"Total soil ice content [mm] = "<< std::setprecision(12) << total_ice_content_sim*1000.0 <<"\n";
  std::cout<<"Soil ice fraction      [-]  = "<< std::setprecision(12) << soil_ice_fraction_sim <<"\n";
    
  passed = test_status > 0 ? "Yes" : "No";
  std::cout<<BLUE<<"\n";
  std::cout<<"\n*********************************************************\n";
  std::cout<<"*************** Summary of the Unit Test ***************\n";
  std::cout<<"*********************************************************\n";
  std::cout<<"Soil ice content error [mm] = "<< err_ice_content_mm <<"\n";
  std::cout<<"Soil ice fraction error [-] = "<< err_ice_fraction <<"\n";
  std::cout<<"Test passed = "<< passed <<"\n";
  std::cout<<RESET<<"\n";

  std::cout<<"\n*********************************************************\n";
  std::cout<<"*********** Testing Calibratable parameters .......\n";
  std::cout<<"\n*********************************************************\n";
  
  double smcmax_set, b_set, satpsi_set;
  ground_temp = 280.15;
  
  model_calib.GetValue("smcmax", &smcmax_set);
  model_calib.GetValue("b", &b_set);
  model_calib.GetValue("satpsi", &satpsi_set);
  model_calib.SetValue("ground_temperature", &ground_temp);
  
  std::cout<<"Initial values | smcmax = "<< smcmax_set <<" , b = "<< b_set <<" , satpsi = "<< satpsi_set <<"\n";

  double smcmax_get, b_get, satpsi_get;
  
  for (int n=0; n<2; n++) {
    smcmax_set += 0.01;
    b_set      += 0.05;
    satpsi_set      += 0.02;
    
    std::cout<<"------------------------------------------------------ \n";
    std::cout<<"Setting | smcmax = "<< smcmax_set <<" , b = "<< b_set <<" , satpsi = "<< satpsi_set <<"\n";
    std::cout<<"------------------------------------------------------ \n";
    model_calib.SetValue("smcmax", &smcmax_set);
    model_calib.SetValue("b", &b_set);
    model_calib.SetValue("satpsi", &satpsi_set);

    model_calib.GetValue("smcmax", &smcmax_get);
    model_calib.GetValue("b", &b_get);
    model_calib.GetValue("satpsi", &satpsi_get);
    std::cout<<"Getting | smcmax = "<< smcmax_get <<" , b = "<< b_get <<" , satpsi = "<< satpsi_get <<"\n";
    std::cout<<"------------------------------------------------------ \n";
    model_calib.Update();
  }
  
  return FAILURE;
}
