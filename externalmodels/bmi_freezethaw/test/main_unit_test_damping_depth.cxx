#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "../../../include/bmi.hxx"
#include "../include/bmi_freezethaw.hxx"
#include "../include/freezethaw.hxx"

#define FAILURE 0
#define VERBOSITY 1

int main(int argc, char *argv[])
{
  BmiFreezeThaw model;

  if (argc != 2) {
    printf("Usage: run_bmifrozensoilcxx CONFIGURATION_FILE\n\n");
    printf("Run the frozensoilcxx model through its BMI with a configuration file.\n");
    return FAILURE;
  }

  std::cout<<"\n**************** BEGIN SoilFreezeThaw (Damping Depth) UNIT TEST *******************\n";

  model.Initialize(argv[1]);
    
  
  

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
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

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  // Test get_time_step()
  dt = model.GetTimeStep();
  if (dt == timestep) {
    if (VERBOSITY)
      std::cout<<" timestep: "<< dt<<"\n";
  }
  else
    test_status &= false;

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
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

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  // Test get_current_time()
  time = model.GetCurrentTime();
  if (VERBOSITY)
    std::cout<<" current time: "<< time<<"\n";

  // Test BMI: GET VALUE FUNCTIONS
  std::cout<<"\n\n************** TEST BMI GETTER SETTER FUNCTIONS********************************\n";

  double time_cur;
    
  // Print current time step - function already tested
  time_cur = model.GetCurrentTime();
  std::cout<<"Current time: "<<time_cur<<"\n";
  
  std::cout<<"********** Input variables ***************** \n";
  // Loop through both input and output variables and call get/set_value_*()
  for (int i=0; i<count_in; i++) {
    std::string var_name = names_in[i];
    std::cout<<"Variable name: "<< var_name <<"\n";
    
    double *var = new double[nz];
    double *dest = new double[nz];
    int indices[] = {0,1,2,3};
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_value() at each timestep
    model.GetValue(var_name, &(var[0]));
    std::cout<<" Get value: "<< var[0] <<"\n";
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_value_at_indices()
    model.GetValueAtIndices(var_name, dest, indices, nz);
    std::cout<<" Get value at indices: " << dest[0]<<"\n";
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_value_ptr()
    double *var_ptr = new double[nz];
    var_ptr= (double*) model.GetValuePtr(var_name);
    std::cout<<" Get value ptr: "<<*var_ptr<<"\n";
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    if (var_name == "soil__moisture_content_total") {
      if (var[0] == soil_MCT[0] && var[1] == soil_MCT[1] && var[2] == soil_MCT[2] && var[3] == soil_MCT[3])
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::stringstream errMsg;
	errMsg << "Soil temperatures should be: "<<soil_T[0]<<" "<<soil_T[1]<<" "<<soil_T[2]<<" "<<soil_T[3]<<"\n";
	throw std::runtime_error(errMsg.str());
      }
    }
    
    // Go ahead and test set_value_*() for last time step here
    //std::cout<<"**Last Time step ** \n";
    // Test BMI set_value()
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    double var_new[] = {0.1,0.2,0.3,0.4};
    double *var_new_up = new double[nz];
    model.SetValue(var_name, &(var_new[0]));
    
    std::cout<<" Set value: "<<var_new[0] <<"\n";
    // get_value to see if changed
    model.GetValue(var_name, &var_new_up[0]);
    std::cout<<" Get value: "<< var_new_up[0] <<"\n";
    
    if (var_name == "soil__moisture_content_liquid") {
      if (var_new[0] == var_new_up[0] && var_new[1] == var_new_up[1] && var_new[2] == var_new_up[2] && var_new[3] == var_new_up[3])
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::stringstream errMsg;
	errMsg << "Getter/Setters are not working properly.\n";
	throw std::runtime_error(errMsg.str());
      }
    }
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test BMI set_value_at_indices()
    double dest_new[] = {0.15,0.25,0.35,0.45};
    double *dest_new_up = new double[nz];
    
    //      int indices[] = {0,1,2,3}; 
    model.SetValueAtIndices(var_name, &indices[0], nz, &dest_new[0]);
    
    std::cout<<" Set value at indices: "<<dest_new[0]<<"\n";
    // get_value_at_indices to see if changed
    model.GetValueAtIndices(var_name, dest_new_up,  &indices[0], nz);
    std::cout<<" Get value at indices: "<<dest_new_up[0]<<"\n";
    if (dest_new[3] == dest_new_up[3])
      test_status &= true;
    else
      test_status &= false;

    // Reset to initial values
    if (var_name == "soil__moisture_content_total") 
      model.SetValue(var_name, &(soil_MCT[0]));
    if (var_name == "soil__moisture_content_liquid") 
      model.SetValue(var_name, &(soil_MCL[0]));
  }
  
  std::cout<<"************* Output variables ***************** \n";
  
  for (int i=0; i<count_out; i++) {
    
    std::string var_name = names_out[i];
    std::cout<<"variable name: "<< var_name <<" "<<test_status<<"\n";
    int len = 0;
    int *indices = NULL;
    
    if (var_name.compare("soil__num_cells") != 0) {
      if (var_name.compare("soil__ice_fraction") == 0 ) {
	len = 1;
	indices = new int[len];
	indices[0] = 0;
      }
      else {
	len = nz;
	indices = new int[len];
	for (int k1=0;k1<len;k1++)
	  indices[k1] = k1;
      }
      
      double *var = new double[len];
      double *dest = new double[len];

      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      // Test get_value() at initial timestep
      model.GetValue(var_name, &(var[0]));
      std::cout<<" Get value: "<< var[0] <<"\n";
      
      if (var_name == "soil__temperature") {
	if (var[0] == soil_T[0] && var[1] == soil_T[1] && var[2] == soil_T[2] && var[3] == soil_T[3])
	  test_status &= true;
	else {
	  test_status &= false;
	  std::string passed = test_status == true ? "Yes" : "No";
	  std::cout<<"Test passed: "<<passed<<"\n";
	  std::stringstream errMsg;
	  errMsg << "Soil temperatures should be: "<<soil_T[0]<<" "<<soil_T[1]<<" "<<soil_T[2]<<" "<<soil_T[3]
		 <<" but are: "<<var[0]<<" "<<var[1]<<" "<<var[2]<<" "<<var[3]<<"\n";
	  throw std::runtime_error(errMsg.str());
	}
      }
      
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      // Test get_value_at_indices()
      model.GetValueAtIndices(var_name, dest, &indices[0], len);
      std::cout<<" Get value at indices: " << dest[0]<<"\n";
      
      if (var_name == "soil__moisture_content_total") {
	if (dest[0] == soil_MCT[0] && dest[1] == soil_MCT[1] && dest[2] == soil_MCT[2] && dest[3] == soil_MCT[3])
	  test_status &= true;
	else {
	  test_status &= false;
	  std::string passed = test_status == true ? "Yes" : "No";
	  std::cout<<"Test passed: "<<passed<<"\n";
	  std::stringstream errMsg;
	  errMsg << "Soil moisture content should be: "<<soil_MCT[0]<<" "<<soil_MCT[1]<<" "<<soil_MCT[2]<<" "<<soil_MCT[3]<<" but are: "<<dest[0]<<" "<<dest[1]<<" "<<dest[2]<<" "<<dest[3]<<"\n";
	  throw std::runtime_error(errMsg.str());
	}
      }
      
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      // Test get_value_ptr()
      double *var_ptr = new double[len];
      var_ptr= (double*) model.GetValuePtr(var_name);
      std::cout<<" Get value ptr: "<<*var_ptr<<"\n";

      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      // Test BMI set_value() 
      double *var_new = new double[len];
      double *var_new_up = new double[len];
      if (var_name == "soil__temperature") {
	for (int i1=0; i1<nz;i1++)
	  var_new[i] = 290+1.5*(i+1);
      }
      else if (var_name == "soil__ice_fraction")
	var_new[0] = 0.05;
      else
	var_new[0] = 0.324;
      model.SetValue(var_name, &(var_new[0]));
	
      std::cout<<" Set value: "<<var_new[0]<<"\n";
      // get_value to see if changed
      model.GetValue(var_name, &var_new_up[0]);
      std::cout<<" Get value: "<< var_new_up[0] <<"\n";

      if (var_name == "soil__temperature") {
	if (var_new[0] == var_new_up[0] && var_new[1] == var_new_up[1] && var_new[2] == var_new_up[2] && var_new[3] == var_new_up[3])
	  test_status &= true;
	else {
	  test_status &= false;
	  std::string passed = test_status == true ? "Yes" : "No";
	  std::cout<<"Test passed: "<<passed<<"\n";
	  std::stringstream errMsg;
	  errMsg << "Getter/Setters are not working properly for output variables.\n";
	  throw std::runtime_error(errMsg.str());
	}
      }

      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      // Test BMI set_value_at_indices()
      double dest_new[] = {0.1,0.21,0.33,0.4};
      if (var_name == "soil__temperature") {
	for (int i1=0; i1<nz;i1++)
	  dest_new[i] = 290+1.5*(i+1);
      }
      
      double *dest_new_up = new double[len];
	
      int indices[] = {0,1,2,3}; 
      model.SetValueAtIndices(var_name, &indices[0], len, &dest_new[0]);
      
      std::cout<<" Set value at indices: "<<dest_new[0]<<"\n";
      // get_value_at_indices to see if changed
      model.GetValueAtIndices(var_name, dest_new_up,  &indices[0], len);
      std::cout<<" Get value at indices: "<<dest_new_up[0]<<"\n";
      int c = 2;
      if (var_name == "soil__ice_fraction") // frozen fraction is just a single number
	c=0;
      if (dest_new[c] == dest_new_up[c])
	test_status &= true;
      else
	test_status &= false;
    }
    else if (var_name.compare("soil__num_cells") == 0) {
      len = 1;
      indices = new int[len];
      indices[0] = 0;
      
      int *var = new int[len];
      
      // Test get_value() at each timestep
      model.GetValue(var_name, &(var[0]));
      std::cout<<" Get value: "<< var[0] <<"\n";
      if (var[0] == nz)
	test_status &= true;
      else
	test_status &= false;
    }
  }
  
  std::string passed = test_status > 0 ? "Yes" : "No";
  std::cout<<"Test passed at time t=0: "<<passed<<"\n";


  int test_nstep=3333;
  std::cout<<"*********** Timestepping now....... (total timesteps): "<< test_nstep<<" **********\n";

  for (int n=0;n<test_nstep;n++) {
    std::cout<<"********* Advancing timestep: "<< (n+1) <<"\n";

    // Test BMI: CONTROL FUNCTION update()
    model_cyc.Update();    
  }

  double solution[] = {272.519, 273.15,274.181,275.154};
  double ice_fraction = 0.0374381;
  double *var_st = new double[nz];
  double *ice_frac = new double[1];
  model_cyc.GetValue("soil__temperature",&var_st[0]);
  model_cyc.GetValue("soil__ice_fraction",&ice_frac[0]);
  double error = 0;
  for (int i1=0; i1<nz;i1++) {
    error += std::pow(var_st[i1] - solution[i1],2.);
    //std::cout<<"SS: "<<var_st[i1]<<" "<<ice_frac[0]<<" "<<*ice_frac<<"\n";
  }
  error = std::pow(error,0.5);
  if (error < 1.e-3)
    test_status &= true;
  else
    test_status &= false;

  double err_frozen_frac_mm = fabs(ice_fraction - ice_frac[0]) *1000; 
  if (err_frozen_frac_mm < 1.e-3)
    test_status &= true;
  else
    test_status &= false;
    
  passed = test_status > 0 ? "Yes" : "No";
  std::cout<<"\n\n*********************************************************\n";
  std::cout<<"*************** Summary of the Unit Test ***************\n";
  std::cout<<"*********************************************************\n";
  std::cout<<"Test passed = "<<passed<<" \nError (L2-norm) =  "<<error<<" \nFrozen fraction = "<<err_frozen_frac_mm<<"\n";
  return FAILURE;
}
