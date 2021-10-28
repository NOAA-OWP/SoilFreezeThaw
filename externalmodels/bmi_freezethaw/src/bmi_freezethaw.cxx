#ifndef BMI_FS_C_INCLUDED
#define BMI_FS_C_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "../../../include/bmi.hxx"
#include "../include/bmi_freezethaw.hxx"
#include "../include/freezethaw.hxx"

void BmiFreezeThaw::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 )
    this->_model = freezethaw::FreezeThaw(config_file);
}


void BmiFreezeThaw::
Update()
{
  this->_model.advance_in_time();
}


void BmiFreezeThaw::
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
    this->_model.dt = frac * dt;
    this->_model.advance_in_time();
    this->_model.dt = dt;
  }
}


void BmiFreezeThaw::
Finalize()
{
  this->_model.~FreezeThaw();
}


int BmiFreezeThaw::
GetVarGrid(std::string name)
{
  if (name.compare("soil__temperature") == 0 || name.compare("soil__mctotal") == 0)
    return 0;
  else if (name.compare("soil__mcliq") == 0 || name.compare("soil__mcice") == 0)
    return 0;
  else
    return -1;
}


std::string BmiFreezeThaw::
GetVarType(std::string name)
{
  if (name.compare("soil__temperature") == 0 || name.compare("soil__mctotal") == 0)
    return "double";
  else if (name.compare("soil__mcliq") == 0 || name.compare("soil__mcice") == 0)
    return "double";
  else
    return "";
}


int BmiFreezeThaw::
GetVarItemsize(std::string name)
{
  if (name.compare("soil__temperature") == 0 || name.compare("soil__mctotal") == 0)
    return sizeof(double);
  else if (name.compare("soil__mcliq") == 0 || name.compare("soil__mcice") == 0)
    return sizeof(double);
  else
    return 0;
}


std::string BmiFreezeThaw::
GetVarUnits(std::string name)
{
  if (name.compare("soil__temperature") == 0)
    return "K";
  else
    return "";
}


int BmiFreezeThaw::
GetVarNbytes(std::string name)
{
  int itemsize;
  int gridsize;

  itemsize = this->GetVarItemsize(name);
  gridsize = this->GetGridSize(this->GetVarGrid(name));
  return itemsize * gridsize;
}


std::string BmiFreezeThaw::
GetVarLocation(std::string name)
{
  if (name.compare("soil__temperature") == 0)
    return "node";
  else if (name.compare("soil__mctotal") == 0)
    return "node";
  else if (name.compare("soil__mcliq") == 0 || name.compare("soil__mcice") == 0)
    return "node";
  else
    return "";
}


void BmiFreezeThaw::
GetGridShape(const int grid, int *shape)
{
  if (grid == 0) {
    shape[0] = this->_model.shape[0];
  }
}


void BmiFreezeThaw::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->_model.spacing[0];
  }
}


void BmiFreezeThaw::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->_model.origin[0];
  }
}


int BmiFreezeThaw::
GetGridRank(const int grid)
{
  if (grid == 0)
    return 1;
  else
    return -1;
}


int BmiFreezeThaw::
GetGridSize(const int grid)
{
  if (grid == 0)
    return this->_model.shape[0];
  else
    return -1;
}


std::string BmiFreezeThaw::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  else
    return "";
}


void BmiFreezeThaw::
GetGridX(const int grid, double *x)
{
  throw NotImplemented();
}


void BmiFreezeThaw::
GetGridY(const int grid, double *y)
{
  throw NotImplemented();
}


void BmiFreezeThaw::
GetGridZ(const int grid, double *z)
{
  throw NotImplemented();
}


int BmiFreezeThaw::
GetGridNodeCount(const int grid)
{
  if (grid == 0)
    return this->_model.shape[0];
  else
    return -1;
}


int BmiFreezeThaw::
GetGridEdgeCount(const int grid)
{
  throw NotImplemented();
}


int BmiFreezeThaw::
GetGridFaceCount(const int grid)
{
  throw NotImplemented();
}


void BmiFreezeThaw::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  throw NotImplemented();
}


void BmiFreezeThaw::
GetGridFaceEdges(const int grid, int *face_edges)
{
  throw NotImplemented();
}


void BmiFreezeThaw::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  throw NotImplemented();
}


void BmiFreezeThaw::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  throw NotImplemented();
}


void BmiFreezeThaw::
GetValue (std::string name, void *dest)
{
  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  nbytes = this->GetVarNbytes(name);
  
  memcpy (dest, src, nbytes);
}


void *BmiFreezeThaw::
GetValuePtr (std::string name)
{
  if (name.compare("soil__temperature") == 0)
    return (void*)this->_model.ST;
  else if (name.compare("soil__mctotal") == 0)
    return (void*)this->_model.SMCT;
  else if (name.compare("soil__mcliq") == 0)
    return (void*)this->_model.SMCLiq;
  else if (name.compare("soil__mcice") == 0)
    return (void*)this->_model.SMCIce;
  else
    return NULL;
}


void BmiFreezeThaw::
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


void BmiFreezeThaw::
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


void BmiFreezeThaw::
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


std::string BmiFreezeThaw::
GetComponentName()
{
  return "The Frozen Soil Model";
}


int BmiFreezeThaw::
GetInputItemCount()
{
  return this->input_var_name_count;
}


int BmiFreezeThaw::
GetOutputItemCount()
{
  return this->output_var_name_count;
}


std::vector<std::string> BmiFreezeThaw::
GetInputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);

  return names;
}


std::vector<std::string> BmiFreezeThaw::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


double BmiFreezeThaw::
GetStartTime () {
  return 0.;
}


double BmiFreezeThaw::
GetEndTime () {
  return this->_model.endtime;
}


double BmiFreezeThaw::
GetCurrentTime () {
  return this->_model.time;
}


std::string BmiFreezeThaw::
GetTimeUnits() {
  return "s";
}


double BmiFreezeThaw::
GetTimeStep () {
  return this->_model.dt;
}

#endif
