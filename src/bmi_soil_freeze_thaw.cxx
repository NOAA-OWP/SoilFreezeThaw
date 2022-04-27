#ifndef BMI_FS_C_INCLUDED
#define BMI_FS_C_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "../bmi/bmi.hxx"
#include "../include/bmi_soil_freeze_thaw.hxx"
#include "../include/soil_freeze_thaw.hxx"
#include <algorithm>

void BmiSoilFreezeThaw::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 )
    this->_model = new soilfreezethaw::SoilFreezeThaw(config_file);
}

void BmiSoilFreezeThaw::
Update()
{
  this->_model->Advance();
}


void BmiSoilFreezeThaw::
UpdateUntil(double t)
{
  double time;
  double dt;

  time = this->GetCurrentTime();
  dt = this->GetTimeStep();

  {
    double n_steps = (t - time) / dt;
    double frac;

    for (int n=0; n<int(n_steps); n++)
      this->Update();

    frac = n_steps - int(n_steps);
    this->_model->dt = frac * dt;
    this->_model->Advance();
    this->_model->dt = dt;
  }
}


void BmiSoilFreezeThaw::
Finalize()
{
  if (this->_model)
    this->_model->~SoilFreezeThaw();
}

int BmiSoilFreezeThaw::
GetVarGrid(std::string name)
{
  if (name.compare("num_cells") == 0 || name.compare("ice_fraction_scheme_bmi") == 0)
    return 0; // int
  else if (name.compare("ground_temperature") == 0 || name.compare("ice_fraction_schaake") == 0 || name.compare("ice_fraction_xinan") == 0)
    return 1; //double
  else if (name.compare("soil_moisture_profile") == 0)
    return 2; // arrays
  else
    return -1;
}


std::string BmiSoilFreezeThaw::
GetVarType(std::string name)
{
  if (name.compare("num_cells") == 0 || name.compare("ice_fraction_scheme_bmi") == 0)
    return "int";
  else  if (name.compare("ground_temperature") == 0 || name.compare("ice_fraction_schaake") == 0 || name.compare("ice_fraction_xinan") == 0)
    return "double";
  else if (name.compare("soil_moisture_profile") == 0)
    return "double";
  else
    return "";
}


int BmiSoilFreezeThaw::
GetVarItemsize(std::string name)
{
  if (name.compare("num_cells") == 0 || name.compare("ice_fraction_scheme_bmi") == 0)
    return sizeof(int);
  else if (name.compare("ground_temperature") == 0 || name.compare("ice_fraction_schaake") == 0 || name.compare("ice_fraction_xinan") == 0)
    return sizeof(double);
  else if (name.compare("soil_moisture_profile") == 0)
    return sizeof(double);
  else
    return 0;
}


std::string BmiSoilFreezeThaw::
GetVarUnits(std::string name)
{
  if (name.compare("ground_temperature") == 0)
    return "K";
  else if (name.compare("ice_fraction_schaake") == 0)
    return "m";
  else
    return "";
}


int BmiSoilFreezeThaw::
GetVarNbytes(std::string name)
{
  int itemsize;
  int gridsize;

  itemsize = this->GetVarItemsize(name);
  gridsize = this->GetGridSize(this->GetVarGrid(name));

  return itemsize * gridsize;
}


std::string BmiSoilFreezeThaw::
GetVarLocation(std::string name)
{
  if (name.compare("ground_temperature") == 0)
    return "node";
  else if (name.compare("ice_fraction_xinan") == 0)
    return "node";
  else if (name.compare("ice_fraction_schaake") == 0 ||  name.compare("num_cells") == 0)
    return "node";
  else if (name.compare("soil_moisture_profile") == 0)
    return "node";
  else
    return "";
}


void BmiSoilFreezeThaw::
GetGridShape(const int grid, int *shape)
{
  if (grid == 0) {
    shape[0] = this->_model->shape[0];
  }
}


void BmiSoilFreezeThaw::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->_model->spacing[0];
  }
}


void BmiSoilFreezeThaw::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->_model->origin[0];
  }
}


int BmiSoilFreezeThaw::
GetGridRank(const int grid)
{
  if (grid == 0)
    return 1;
  else
    return -1;
}


int BmiSoilFreezeThaw::
GetGridSize(const int grid)
{
  if (grid == 2) // for arrays
    return this->_model->shape[0];
  if (grid == 0 || grid == 1) // for scalars
    return 1;
  else
    return -1;
}


std::string BmiSoilFreezeThaw::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  else
    return "";
}


void BmiSoilFreezeThaw::
GetGridX(const int grid, double *x)
{
  throw NotImplemented();
}


void BmiSoilFreezeThaw::
GetGridY(const int grid, double *y)
{
  throw NotImplemented();
}


void BmiSoilFreezeThaw::
GetGridZ(const int grid, double *z)
{
  throw NotImplemented();
}


int BmiSoilFreezeThaw::
GetGridNodeCount(const int grid)
{
  if (grid == 0)
    return this->_model->shape[0];
  else
    return -1;
}


int BmiSoilFreezeThaw::
GetGridEdgeCount(const int grid)
{
  throw NotImplemented();
}


int BmiSoilFreezeThaw::
GetGridFaceCount(const int grid)
{
  throw NotImplemented();
}


void BmiSoilFreezeThaw::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  throw NotImplemented();
}


void BmiSoilFreezeThaw::
GetGridFaceEdges(const int grid, int *face_edges)
{
  throw NotImplemented();
}


void BmiSoilFreezeThaw::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  throw NotImplemented();
}


void BmiSoilFreezeThaw::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  throw NotImplemented();
}


void BmiSoilFreezeThaw::
GetValue (std::string name, void *dest)
{
  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  nbytes = this->GetVarNbytes(name);
  memcpy (dest, src, nbytes);
}


void *BmiSoilFreezeThaw::
GetValuePtr (std::string name)
{
  /*
  if (name.compare("temperature") == 0)
    return (void*)this->_model->ST;
  else if (name.compare("moisture_content_total") == 0)
    return (void*)this->_model->SMCT;
  else if (name.compare("moisture_content_liquid") == 0)
    return (void*)this->_model->SMCLiq;
  else if (name.compare("moisture_content_ice") == 0)
    return (void*)this->_model->SMCIce;
  else 
  */
  if (name.compare("soil_moisture_profile") == 0)
    return (void*)this->_model->soil_moisture_content;
  else if (name.compare("ground_temperature") == 0 )
    return (void*)(&this->_model->ground_temp);
  else if (name.compare("num_cells") == 0)
    return (void*)(&this->_model->ncells);
  else if (name.compare("ice_fraction_schaake") == 0) {
    return (void*)(&this->_model->ice_fraction_schaake);
  }
  else if (name.compare("ice_fraction_xinan") == 0){
    int val = this->_model->ice_fraction_xinan;
    return (void*)(&this->_model->ice_fraction_xinan);
  }
  else if (name.compare("ice_fraction_scheme_bmi") == 0)
    return (void*)(&this->_model->ice_fraction_scheme_bmi);
  else {
    std::stringstream errMsg;
    errMsg << "variable "<< name << " does not exist";
    throw std::runtime_error(errMsg.str());
    return NULL;
  }
}


void BmiSoilFreezeThaw::
GetValueAtIndices (std::string name, void *dest, int *inds, int len)
{
  void * src = NULL;

  src = this->GetValuePtr(name);

  if (src) {
    int i;
    int itemsize = 0;
    int offset;
    char *ptr;

    itemsize = this->GetVarItemsize(name);

    for (i=0, ptr=(char *)dest; i<len; i++, ptr+=itemsize) {
      offset = inds[i] * itemsize;
      memcpy(ptr, (char *)src + offset, itemsize);
    }
  }
}


void BmiSoilFreezeThaw::
SetValue (std::string name, void *src)
{
  void * dest = NULL;
  
  dest = this->GetValuePtr(name);
  
  if (dest) {
    int nbytes = 0;
    nbytes = this->GetVarNbytes(name);
    memcpy(dest, src, nbytes);
  }

}


void BmiSoilFreezeThaw::
SetValueAtIndices (std::string name, int * inds, int len, void *src)
{
  void * dest = NULL;

  dest = this->GetValuePtr(name);

  if (dest) {
    int i;
    int itemsize = 0;
    int offset;
    char *ptr;

    itemsize = this->GetVarItemsize(name);

    for (i=0, ptr=(char *)src; i<len; i++, ptr+=itemsize) {
      offset = inds[i] * itemsize;
      memcpy((char *)dest + offset, ptr, itemsize);
    }
  }
}


std::string BmiSoilFreezeThaw::
GetComponentName()
{
  return "Soil Freeze Thaw Model";
}


int BmiSoilFreezeThaw::
GetInputItemCount()
{
  return this->input_var_name_count;
}


int BmiSoilFreezeThaw::
GetOutputItemCount()
{
  return this->output_var_name_count;
}


std::vector<std::string> BmiSoilFreezeThaw::
GetInputVarNames()
{
  std::vector<std::string> names;

  std::vector<std::string>* names_m = _model->InputVarNamesModel();
  
  for (int i=0; i<this->input_var_name_count; i++) {
    if (std::find(names_m->begin(), names_m->end(), this->input_var_names[i]) != names_m->end()) {
      names.push_back(this->input_var_names[i]);
    }
  }
  
  return names;
}


std::vector<std::string> BmiSoilFreezeThaw::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->output_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


double BmiSoilFreezeThaw::
GetStartTime () {
  return 0.;
}


double BmiSoilFreezeThaw::
GetEndTime () {
  return this->_model->endtime;
}


double BmiSoilFreezeThaw::
GetCurrentTime () {
  return this->_model->time;
}


std::string BmiSoilFreezeThaw::
GetTimeUnits() {
  return "s";
}


double BmiSoilFreezeThaw::
GetTimeStep () {
  return this->_model->dt;
}

#endif
