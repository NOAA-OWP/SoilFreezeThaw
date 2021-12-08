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
  BmiFreezeThaw model,model_cyc;

  if (argc != 2) {
    printf("Usage: run_bmifrozensoilcxx CONFIGURATION_FILE\n\n");
    printf("Run the frozensoilcxx model through its BMI with a configuration file.\n");
    return FAILURE;
  }

  std::cout<<"\n**************** BEGIN SoilFreezeThaw BMI UNIT TEST *******************\n";

  model.Initialize(argv[1]);
  model_cyc.Initialize(argv[1]);

  std::cout<<"\n**************** TEST VALUES ************************************\n";
  int nz = 4;
  double endtime = 1035.91*86400.;
  double timestep = 3600;
  bool test_status = true;
  int num_input_vars = 2;
  int num_output_vars = 6;
  int nbytes_input = nz * sizeof(double);
  double soil_MCT[] = {0.36,0.39,0.41,0.43}; // total_moisture_content
  double soil_MCL[] = {0.36,0.39,0.41,0.43}; //liquid_moisture_content
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
    std::cout<<"Why? Number of input variables are different."<<"\n";
    abort();
  }

  // Test GetInputVarNames 
  names_in = model.GetInputVarNames();
  if (VERBOSITY) {
    std::cout<<"Input variable names \n";
    for (int i=0; i<count_in; i++)
      std::cout<<names_in[i]<<"\n";
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
      std::cout<<names_out[i]<<"\n";
  }
  if (count_out == num_output_vars)
    test_status &= true;
  else {
    test_status &= false;
    std::string passed = test_status == true ? "Yes" : "No";
    std::cout<<"Test passed: "<<passed<<"\n";
    std::cout<<"Why? Number of output variables are different."<<"\n";
    abort();
  }
  
  // Test BMI: VARIABLE INFORMATION FUNCTIONS
  std::cout<<"\n**************** TEST BMI VARIABLE INFORMATION FUNCTIONS\n*****************************************\n";
  int grid, itemsize, nbytes;
  std::string location;
  std::string units;
  std::string vartype;
  
  // Loop through both input and output variables and call GetVar* functions
  for (int i=0; i<count_in; i++) {
    std::string var_name = names_in[i];
    if (VERBOSITY)
      std::cout<<"Input var_name: "<< var_name <<"\n";

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_grid()
    grid = model.GetVarGrid(var_name);
    if (VERBOSITY)
      std::cout<<"Grid: "<< grid <<"\n";

    if (grid >=0)
      test_status &= true;
    else
      test_status &= false;

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_itemsize()
    itemsize = model.GetVarItemsize(var_name);
    if (VERBOSITY)
      std::cout<<"Itemsize: "<< itemsize <<"\n";
    if (itemsize >0)
      test_status &= true;
    else
      test_status &= false;

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_location()
    location = model.GetVarLocation(var_name);
    if ( location == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" location: "<< location<<"\n";
    if (location == "")
      test_status &= false;
    else
      test_status &= true;
	
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_units()
    units = model.GetVarUnits(var_name);
    if (VERBOSITY)
      std::cout<<" units: ["<< units <<"]\n";

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_type()
    vartype = model.GetVarType(var_name);
    if (VERBOSITY)
      std::cout<<" type: "<< vartype <<"\n";
    if (location == "")
      test_status &= false;
    else
      test_status &= true;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // get_var_nbytes()
    nbytes = model.GetVarNbytes(var_name);
    if (nbytes == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<" nbytes: "<< nbytes<<"\n";

    if (var_name == "soil__moisture_content_total") {
      if (nbytes == nbytes_input)
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::cout<<"Why? number of bytes for input var"<<var_name<< " should be "<<nbytes_input<<"\n";
	abort();
      }
    }
  }

  if (VERBOSITY)
    std::cout<<"\n*****************************************\n";
  
  for (int i=0; i<count_out; i++) {
    std::string var_name = names_out[i];
    if (VERBOSITY)
      std::cout<<"Output var_name: "<< var_name <<"\n";
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_grid() 
    grid = model.GetVarGrid(var_name);
    if (grid == -1) return -1;
    if (VERBOSITY)
      std::cout<<"Grid: "<< grid <<"\n";

    if (grid >=0)
      test_status &= true;
    else
      test_status &= false;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_itemsize()
    itemsize = model.GetVarItemsize(var_name);
    if (itemsize == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<"Itemsize: "<< itemsize <<"\n";

    if (itemsize >0)
      test_status &= true;
    else
      test_status &= false;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_location()
    location = model.GetVarLocation(var_name);
    if ( location == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" location:"<< location<<"\n";

    if (location == "")
      test_status &= false;
    else
      test_status &= true;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_units()
    units = model.GetVarUnits(var_name);
    if (VERBOSITY)
      std::cout<<" units: ["<< units <<"]\n";
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_type()
    vartype = model.GetVarType(var_name);
    if (vartype == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" type: "<< vartype <<"\n";
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // get_var_nbytes()
    nbytes = model.GetVarNbytes(var_name);
    if (nbytes == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<" nbytes: "<< nbytes<<"\n";

    if (var_name == "soil__temperature") {
      if (nbytes == nbytes_input)
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::cout<<"Why? number of bytes for out var"<<var_name<< " should be "<<nbytes_input<<"\n";
	abort();
      }
    }
    
  }

  // Test BMI: MODEL GRID FUNCTIONS
  std::cout<<"\n \n**************** TEST BMI GRID FUNCTIONS***********************\n";
  int grid_id = 0; //TODO: impliment for multiple grids, for now we know 0
  int grid_rank, grid_size;
  std::string grid_type;

  if (VERBOSITY)  
    std::cout<<" Grid id "<< grid_id<<"\n";

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  // Test get_grid_rank()
  grid_rank = model.GetGridRank(grid_id);
  if (grid_rank == FAILURE) return FAILURE;
  if (VERBOSITY)
    std::cout<<" rank: "<<grid_rank<<"\n";

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  // Test get_grid_size
  grid_size = model.GetGridSize(grid_id);
  if (grid_size == nz) {
    test_status &= true;
    if (VERBOSITY)
      std::cout<<"grid size: "<<grid_size<<"\n";
  }
  else {
    test_status &= false;
    std::string passed = test_status == true ? "Yes" : "No";
    std::cout<<"Test passed: "<<passed<<"\n";
    std::cout<<"Why? Grid size of should be "<<nz<<"\n";
    abort();
  }

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*
  // Test get_grid_type
  grid_type = model.GetGridType(grid_id);
  if (VERBOSITY)
  std::cout<<" grid type: "<<grid_type<<"\n";  
  std::cout<<"Test status: "<<test_status<<"\n";
  */

  // Test BMI: TIME FUNCTIONS
  std::cout<<"\n \n**************** TEST BMI TIME FUNCTIONS***********************\n";
  double time = 0.0;
  double dt = 0.0;
  std::string units_time;

  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  // Test get_start_time()
  time = model.GetStartTime();
  if (VERBOSITY)
    std::cout<<" start time: "<< time<<"\n";

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
    std::cout<<"Why? End time should be ["<<endtime/3600.<<" hours] \n";
    abort();
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
    std::cout<<"Why? Time units should be seconds [s], but the model returned"<<units_time<<"\n";
    abort();
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
	std::cout<<"Why? Soil temperatures should be: "<<soil_T[0]<<" "<<soil_T[1]<<" "<<soil_T[2]<<" "<<soil_T[3]<<"\n";
	abort();
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
	std::cout<<"Why? Getter/Setters are not working properly.\n";
	abort();
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
	  std::cout<<"Why? Soil temperatures should be: "<<soil_T[0]<<" "<<soil_T[1]<<" "<<soil_T[2]<<" "<<soil_T[3]<<" ";
	  std::cout<<"but it is: "<<var[0]<<" "<<var[1]<<" "<<var[2]<<" "<<var[3]<<"\n";
	  abort();
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
	  std::cout<<"Why? Soil moisture content should be: "<<soil_MCT[0]<<" "<<soil_MCT[1]<<" "<<soil_MCT[2]<<" "<<soil_MCT[3]<<" ";
	  std::cout<<"but it is: "<<dest[0]<<" "<<dest[1]<<" "<<dest[2]<<" "<<dest[3]<<"\n";
	  abort();
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
	  std::cout<<"Why? Getter/Setters are not working properly for output variables.\n";
	  abort();
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
  std::cout<<"*************** Summary of the Unit Test ***************\n";
  std::cout<<"Test passed = "<<passed<<" \nError (L2-norm) =  "<<error<<" \nFrozen fraction = "<<err_frozen_frac_mm<<"\n";
  return FAILURE;
}
