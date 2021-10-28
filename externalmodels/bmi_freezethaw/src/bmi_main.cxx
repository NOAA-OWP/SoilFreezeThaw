#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "../../../include/bmi.hxx"
#include "../include/bmi_freezethaw.hxx"
#include "../include/freezethaw.hxx"


int main(int argc, char *argv[])
{
  BmiFreezeThaw model;

  if (argc != 2) {
    printf("Usage: run_bmifrozensoilcxx CONFIGURATION_FILE\n\n");
    printf("Run the frozensoilcxx model through its BMI with a configuration file.\n");
    printf("Output is written to the file `bmifrozensoilcxx.out`.\n");
    return bmi::BMI_SUCCESS;
  }
  
  FILE *fp = fopen("bmifrozensoilcxx.out", "w");
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
    std::string var_name = "soil__temperature";
    std::string var_name1 = "soil__mcliq";
    std::string var_name2 = "soil__mcice";
    int grid, rank, *shape;
    double *var = NULL;
    double *var1 = NULL;
    double *var2 = NULL;
    double time = 0.0;
    double end_time;

    fprintf(fp, "variable = %s\n", var_name.c_str());
    fprintf(fp, "variable = %s\n", var_name1.c_str());
    fprintf(fp, "variable = %s\n", var_name2.c_str());
    grid = model.GetVarGrid(var_name);

    rank = model.GetGridRank(grid);
    fprintf(fp, "rank = %d\n", rank);
    shape = new int[rank];
    model.GetGridShape(grid, shape);

    fprintf(fp, "shape = %d x %d x %d\n", shape[0],1,1);

    var = (double *)model.GetValuePtr(var_name);
    var1 = (double *)model.GetValuePtr(var_name1);
    var2 = (double *)model.GetValuePtr(var_name2);

    end_time = model.GetEndTime();
    while (time < end_time) {

      time = model.GetCurrentTime();
      fprintf(fp, "\nTime = %f\n", time);

      for (int i=0; i < shape[0]; i++) {
	fprintf(fp, "%6.4e, %6.4e, %6.4e", var[i],var1[i], var2[i]);
	fprintf(fp, "\n");	
      }
      
      model.Update();
    }
  }

  fprintf(fp, "Finalizing... ");
  model.Finalize();
  fprintf(fp, "done\n");

  fclose(fp);

  return bmi::BMI_SUCCESS;
}
