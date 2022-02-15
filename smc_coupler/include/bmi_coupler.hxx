#ifndef BMI_COUPLER_H_INCLUDED
#define BMI_COUPLER_H_INCLUDED

using namespace std;

#include <string.h>
#include "../../bmi/bmi.hxx"
#include "smc_profile.hxx"

namespace coupler {
class NotImplemented : public std::logic_error {
  public:
  NotImplemented() : std::logic_error("Not Implemented") { };
};

}
class BmiCoupler : public bmixx::Bmi {
  public:
    BmiCoupler() {
      this->input_var_names[0] = "soil__storage";
      this->input_var_names[1] = "soil__storage_change";
      this->input_var_names[2] = "soil__water_table";
      this->input_var_names[3] = "soil__moisture_content_layered";
      
      this->output_var_names[0] = "soil__moisture_content_total";
    };

  void Initialize(std::string config_file);
    
  void Update();
  void UpdateUntil(double time);
    void Finalize();

    std::string GetComponentName();
    int GetInputItemCount();
    int GetOutputItemCount();
    std::vector<std::string> GetInputVarNames();
    std::vector<std::string> GetOutputVarNames();

    int GetVarGrid(std::string name);
    std::string GetVarType(std::string name);
    int GetVarItemsize(std::string name);
    std::string GetVarUnits(std::string name);
    int GetVarNbytes(std::string name);
    std::string GetVarLocation(std::string name);

    double GetCurrentTime();
    double GetStartTime();
    double GetEndTime();
    std::string GetTimeUnits();
    double GetTimeStep();

    void GetValue(std::string name, void *dest);
    void *GetValuePtr(std::string name);
    void GetValueAtIndices(std::string name, void *dest, int *inds, int count);

    void SetValue(std::string name, void *src);
    void SetValueAtIndices(std::string name, int *inds, int len, void *src);

    int GetGridRank(const int grid);
    int GetGridSize(const int grid);
    std::string GetGridType(const int grid);

    void GetGridShape(const int grid, int *shape);
    void GetGridSpacing(const int grid, double *spacing);
    void GetGridOrigin(const int grid, double *origin);

    void GetGridX(const int grid, double *x);
    void GetGridY(const int grid, double *y);
    void GetGridZ(const int grid, double *z);

    int GetGridNodeCount(const int grid);
    int GetGridEdgeCount(const int grid);
    int GetGridFaceCount(const int grid);

    void GetGridEdgeNodes(const int grid, int *edge_nodes);
    void GetGridFaceEdges(const int grid, int *face_edges);
    void GetGridFaceNodes(const int grid, int *face_nodes);
    void GetGridNodesPerFace(const int grid, int *nodes_per_face);
  private:
    smc_profile::SMCProfile _model;
    static const int input_var_name_count = 4;
    static const int output_var_name_count = 1;

    std::string input_var_names[4];
    std::string output_var_names[1];
};

#endif
