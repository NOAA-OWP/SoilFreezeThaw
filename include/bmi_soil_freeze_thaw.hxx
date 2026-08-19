#ifndef BMI_SFT_H_INCLUDED
#define BMI_SFT_H_INCLUDED

using namespace std;

#include <string.h>
#include "../bmi/bmi.hxx"
#include "soil_freeze_thaw.hxx"
#include "vecbuf.hpp"

#include <boost/serialization/serialization.hpp>

class NotImplemented : public std::logic_error {
  public:
  NotImplemented() : std::logic_error("Not Implemented in SFT") { };
};


class BmiSoilFreezeThaw : public bmi::Bmi {
  public:
    BmiSoilFreezeThaw();
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
    friend class boost::serialization::access;
    soilfreezethaw::SoilFreezeThaw* state;
    vecbuf<char> m_serialized_vec;
    uint64_t m_serialized_length; // can theoretically always be derived from the vec's size, but needed for having a stable location for GetValuePtr
    static const int input_var_name_count  = 2;
    static const int output_var_name_count = 6;
    static const int calib_var_name_count  = 3;

    std::string input_var_names[input_var_name_count];
    std::string output_var_names[output_var_name_count];
    std::string calib_var_names[calib_var_name_count];
    std::string verbosity;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version);
    void new_serialized();
    void load_serialized(char* data);
    void clear_serialized();
};

#ifdef NGEN
extern "C"
{

    /**
    * Construct this BMI instance as a normal C++ object, to be returned to the framework.
    *
    * @return A pointer to the newly allocated instance.
    */
  BmiSoilFreezeThaw *bmi_model_create()
  {
    return new BmiSoilFreezeThaw();
  }
  
    /**
     * @brief Destroy/free an instance created with @see bmi_model_create
     * 
     * @param ptr 
     */
  void bmi_model_destroy(BmiSoilFreezeThaw *ptr)
  {
    delete ptr;
  }

}
#endif

#endif
