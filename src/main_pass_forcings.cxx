#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/cfe.h"
#include "../include/bmi.h"
#include "../include/bmi_cfe.h"

#include "../forcing_code/include/aorc.h"
#include "../forcing_code/include/bmi_aorc.h"

#include "../include/bmi.hxx"
#include "../externalmodels/bmi_freezethawcxx/include/bmi_freezethaw.hxx"
#include "../externalmodels/bmi_freezethawcxx/include/freezethaw.hxx"


/***************************************************************
    Function to pass the ice fraction from Freeze-thaw model to CFE using BMI.
***************************************************************/
void pass_icefraction_from_ftm_to_cfe(Bmi *cfe_bmi_model, BmiFreezeThaw ftm_bmi_model){
  
  /********************************************************************
        TODO: Get variable names through BMI, then loop through those
              so we don't re-write the get/set functions over and over
  ********************************************************************/
  int *nz_ = new int[1];
  double *ice_fraction_v = new double[1];
  
  ftm_bmi_model.GetValue("soil__num_layers", &(nz_[0]));
  ftm_bmi_model.GetValue("soil__moisture_content_ice_bulk", &(ice_fraction_v[0]));
  int nz = *nz_;
  double ice_fraction = *ice_fraction_v;
  /*
  double *smct_b = new double[nz];
  ftm_bmi_model.GetValue("soil__moisture_content_total_bulk", &(smct_b[0]));
  std::cout<<"soil mc liq= "<< *smct_b<<"\n";
  */
  cfe_bmi_model->set_value(cfe_bmi_model, "soil__ice_fraction", &(ice_fraction_v[0]));
  
}

/***************************************************************
    Function to pass the ice fraction from Freeze-thaw model to CFE using BMI.
***************************************************************/
void pass_smc_from_cfe_to_ftm(Bmi *cfe_bmi_model, BmiFreezeThaw ftm_bmi_model){
  
  /********************************************************************
        TODO: Get variable names through BMI, then loop through those
              so we don't re-write the get/set functions over and over
  ********************************************************************/
  double time = ftm_bmi_model.GetCurrentTime () ;
  double end_time = ftm_bmi_model.GetEndTime();

  double *smct_v = new double[8];

  cfe_bmi_model->get_value(cfe_bmi_model, "SMCT", &smct_v[0]);
  //  std::cout<<"soil mc from cfe = "<<smct_v[3] <<"\n";
  ftm_bmi_model.SetValue("soil__moisture_content_total", &(smct_v[0]));
  
}

/***************************************************************
    Function to pass the forcing data from AORC to CFE using BMI.
    This requires a lot of getters and setters, 
    so no need to clutter up main program
***************************************************************/
void pass_forcing_from_aorc_to_cfe(Bmi *cfe_bmi_model, Bmi *aorc_bmi_model){

    /********************************************************************
        TODO: Get variable names through BMI, then loop through those
              so we don't re-write the get/set functions over and over
    ********************************************************************/

    double *var = NULL;
    var = (double*) malloc (sizeof (double)*1);
    
    aorc_bmi_model->get_value(aorc_bmi_model, "atmosphere_water__liquid_equivalent_precipitation_rate", &(var[0]));
    
    cfe_bmi_model->set_value(cfe_bmi_model, "atmosphere_water__liquid_equivalent_precipitation_rate", &(var[0]));

}


/************************************************************************
    This main program is a mock framwork.
    This is not part of BMI, but acts as the driver that calls the model.
************************************************************************/
int
 main(int argc, const char *argv[])
{

  /************************************************************************
      A configuration file is required for running this model through BMI
  ************************************************************************/
  if(argc<=1){
    printf("make sure to include a path to the CFE config file\n");
    exit(1);
  }
  if(argc<=2){
    printf("make sure to include a path to the AORC config file\n");
    exit(1);
  }

  /************************************************************************
      allocating memory to store the entire BMI structure for CFE and AORC
  ************************************************************************/
  printf("allocating memory to store entire BMI structure for CFE\n");
  Bmi *cfe_bmi_model = (Bmi *) malloc(sizeof(Bmi));
  printf("allocating memory to store entire BMI AORC structure\n");
  Bmi *aorc_bmi_model = (Bmi *) malloc(sizeof(Bmi));
  
  BmiFreezeThaw ftm_bmi_model;
  
  /************************************************************************
      Registering the BMI model for CFE and AORC
  ************************************************************************/
  printf("Registering BMI CFE model\n");
  register_bmi_cfe(cfe_bmi_model);
  printf("Registering BMI AORC model\n");
  register_bmi_aorc(aorc_bmi_model);

  /************************************************************************
      Initializing the BMI model for CFE and AORC and Freeze-thaw model
  ************************************************************************/
  printf("Initializeing BMI CFE model\n");

  const char *cfg_file = argv[1];
  cfe_bmi_model->initialize(cfe_bmi_model, cfg_file);
  printf("Initializeing BMI AORC model\n");

  const char *cfg_file_aorc = argv[2];
  aorc_bmi_model->initialize(aorc_bmi_model, cfg_file_aorc);
  printf("Initializeing BMI FTM model\n");

  const char *cfg_file_ftm = argv[3];
  ftm_bmi_model.Initialize(cfg_file_ftm);

  /************************************************************************
    Get the information from the configuration here in Main
  ************************************************************************/
  printf("Get the information from the configuration here in Main\n");
  cfe_state_struct *cfe_model_data;
  cfe_model_data = (cfe_state_struct *) cfe_bmi_model->data;
  aorc_model *aorc;
  aorc = (aorc_model *) aorc_bmi_model->data;

  /************************************************************************
    This is the basic process for getting two things to talk through BMI
    1. Update the AORC forcing data
    2. Getting forcing from AORC and setting forcing for CFE
    3. Update the CFE model.
  ************************************************************************/

  /************************************************************************
    Now loop through time and call the models with the intermediate get/set
  ************************************************************************/
  printf("looping through and calling updata\n");
  if (cfe_model_data->verbosity > 0)
    print_cfe_flux_header();
  for (int i = 0; i < 5; i++){
    aorc_bmi_model->update(aorc_bmi_model);                         // Update model 1
    pass_forcing_from_aorc_to_cfe(cfe_bmi_model, aorc_bmi_model);   // Get and Set values
    pass_icefraction_from_ftm_to_cfe(cfe_bmi_model, ftm_bmi_model);
    if (cfe_model_data->aorc.precip_kg_per_m2 != aorc->aorc.precip_kg_per_m2){
      printf("Precip values do not match\n");
      printf("precip value from AORC is %lf\n", aorc->aorc.precip_kg_per_m2);
      printf("precip value from CFE is %lf\n", cfe_model_data->aorc.precip_kg_per_m2);
    }

    cfe_bmi_model->update(cfe_bmi_model);                           // Update model 2
  
    if (cfe_model_data->verbosity > 0)
      print_cfe_flux_at_timestep(cfe_model_data);
    
    pass_smc_from_cfe_to_ftm(cfe_bmi_model, ftm_bmi_model);
    //updated_ftm_smc_soil_temperatue(cfe_bmi_model, ftm_bmi_model);
    ftm_bmi_model.Update();
  }

  // Run the Mass Balance check
  mass_balance_check(cfe_model_data);

  /************************************************************************
    Finalize both the CFE and AORC bmi models
  ************************************************************************/
  printf("Finalizing BFE and AORC models\n");
  cfe_bmi_model->finalize(cfe_bmi_model);
  aorc_bmi_model->finalize(aorc_bmi_model);

  return 0;
}

