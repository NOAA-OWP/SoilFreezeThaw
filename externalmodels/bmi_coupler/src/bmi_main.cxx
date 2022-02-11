#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "../../../include/bmi.hxx"
#include "../include/bmi_coupler.hxx"
#include "../include/smc_profile.hxx"

#define SUCCESS 0
int main(int argc, char *argv[])
{
  BmiCoupler model;
  
  if (argc != 2) {
    printf("Usage: run_bmifrozensoilcxx CONFIGURATION_FILE\n\n");
    printf("Run the frozensoilcxx model through its BMI with a configuration file.\n");
    printf("Output is written to the file `bmifrozensoilcxx.out`.\n");
    return SUCCESS;
  }

  FILE *fp = fopen("bmi_file.out", "w");
  fprintf(fp, "Configuration file = %s\n", argv[1]);
  fprintf(fp, "Initializing... ");
  
  model.Initialize(argv[1]);
  
  fprintf(fp, "done\n");

  {
    std::string model_name;
    model_name = model.GetComponentName();
    fprintf(fp, "%s\n", model_name.c_str());
  }

  {
    std::string var_name_s = "soil__storage";
    std::string var_name_sc = "soil__storage_change";
    std::string var_name_wt = "soil__water_table";
    std::string var_name_smc = "soil__moisture_content_total";
    
    int grid, rank, *shape;
    double *var_s = NULL;
    double *var_sc = NULL;
    double *var_wt = NULL;

    fprintf(fp, "variable = %s\n", var_name_s.c_str());
    fprintf(fp, "variable = %s\n", var_name_sc.c_str());
    fprintf(fp, "variable = %s\n", var_name_wt.c_str());
    fprintf(fp, "variable = %s\n", var_name_smc.c_str());
    
    grid = model.GetVarGrid(var_name_s);

    rank = model.GetGridRank(grid);
    fprintf(fp, "rank = %d\n", rank);
    shape = new int[rank];
    model.GetGridShape(grid, shape);

    fprintf(fp, "shape = %d x %d x %d\n", shape[0],1,1);

    // Set values
    double storage_m = 0.526328;
    double storage_change_m = -0.000472;
    double *storage_m_ptr = &storage_m;
    double *storage_change_m_ptr = &storage_change_m;
    
    model.SetValue(var_name_s,storage_m_ptr);

  
    model.SetValue(var_name_sc,storage_change_m_ptr);

    
    var_s = (double *)model.GetValuePtr(var_name_s);
    var_sc = (double *)model.GetValuePtr(var_name_sc);

    std::cout<<"storage: "<<*var_s<<"\n";
    std::cout<<"storage change: "<<*var_sc<<"\n";

    var_wt = (double *)model.GetValuePtr(var_name_wt);
    
    std::cout<<"water table: "<<*var_wt<<"\n";

    model.Update();
    // unit test
    double SMCT[] ={0.322036, 0.33341, 0.367307, 0.439};
    // Get values
    double *var_smc = new double[4];
    
    model.GetValue(var_name_smc,&var_smc[0]);


    for (int i=0; i < shape[0]; i++) {
      std::cout<<"Main: "<<var_smc[i]<<" "<<SMCT[i]<<" "<<abs(var_smc[i] - SMCT[i])<<"\n";
      assert (abs(var_smc[i] - SMCT[i]) < 1.E-6);     
      fprintf(fp, "%6.4e", var_smc[i]);
      fprintf(fp, "\n");
    }
    
  }
  
  fprintf(fp, "Finalizing... ");

  model.Finalize();
  fprintf(fp, "done\n");
  fclose(fp);
  return SUCCESS;
}
