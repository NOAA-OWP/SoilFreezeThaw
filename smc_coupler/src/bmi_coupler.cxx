#ifndef BMI_COUPLER_C_INCLUDED
#define BMI_COUPLER_C_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "../../bmi/bmi.hxx"
#include "../include/bmi_coupler.hxx"
#include "../include/smc_profile.hxx"

void BmiCoupler::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 )
    this->_model = smc_profile::SMCProfile(config_file);
}


void BmiCoupler::
Update()
{
  this->_model.SMPVertical();
}


void BmiCoupler::
UpdateUntil(double t)
{
  this->_model.SMPVertical();
}


void BmiCoupler::
Finalize()
{
  this->_model.~SMCProfile();
}


int BmiCoupler::
GetVarGrid(std::string name)
{
  if (name.compare("soil__storage") == 0 || name.compare("soil__storage_change") == 0 || name.compare("soil__water_table") == 0)
    return 0;
  else if (name.compare("soil__moisture_content_total") == 0 || name.compare("soil__moisture_content_layered") == 0)
    return 1;
  else if (name.compare("soil__smc_profile_option_bmi") == 0)
    return 2;
  else 
    return -1;
}


std::string BmiCoupler::
GetVarType(std::string name)
{
   if (name.compare("soil__storage") == 0 || name.compare("soil__storage_change") == 0 || name.compare("soil__water_table") == 0)
    return "double";
  else if (name.compare("soil__moisture_content_total") == 0 || name.compare("soil__moisture_content_layered") == 0)
    return "double";
  else if (name.compare("soil__smc_profile_option_bmi") == 0)
    return "int";
  else
    return "";
}


int BmiCoupler::
GetVarItemsize(std::string name)
{
  if (name.compare("soil__storage") == 0 || name.compare("soil__storage_change") == 0 || name.compare("soil__water_table") == 0)
    return sizeof(double);
  else if (name.compare("soil__moisture_content_total") == 0 || name.compare("soil__moisture_content_layered") == 0)
    return sizeof(double);
  else if (name.compare("soil__smc_profile_option_bmi") == 0)
    return sizeof(int);
  else
    return 0;
}


std::string BmiCoupler::
GetVarUnits(std::string name)
{
  if (name.compare("soil__storage") == 0 || name.compare("soil__storage_change") == 0 || name.compare("soil__water_table") == 0)
    return "m";
  else if (name.compare("soil__moisture_content_total") == 0 || name.compare("soil__moisture_content_layered") == 0)
    return "";
  else
    return "";
}


int BmiCoupler::
GetVarNbytes(std::string name)
{
  int itemsize;
  int gridsize;

  itemsize = this->GetVarItemsize(name);
  gridsize = this->GetGridSize(this->GetVarGrid(name));
  return itemsize * gridsize;
}


std::string BmiCoupler::
GetVarLocation(std::string name)
{
  if (name.compare("soil__storage") == 0 || name.compare("soil__storage_change") == 0)
    return "domain";
  else if (name.compare("soil__water_table") == 0)
    return "node";
  else if (name.compare("soil__moisture_content_total") == 0)
    return "domain";
  else if (name.compare("soil__moisture_content_layered") == 0)
    return "layers";
  else
    return "";
}


void BmiCoupler::
GetGridShape(const int grid, int *shape)
{
  if (grid == 0) {
    shape[0] = this->_model.shape[0];
  }
}


void BmiCoupler::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->_model.spacing[0];
  }
}


void BmiCoupler::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->_model.origin[0];
  }
}


int BmiCoupler::
GetGridRank(const int grid)
{
  if (grid == 0)
    return 1;
  else
    return -1;
}


int BmiCoupler::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 2)
    return 1;
  else if (grid == 1)
    return this->_model.shape[0];
  else
    return -1;
}


std::string BmiCoupler::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  else
    return "";
}


void BmiCoupler::
GetGridX(const int grid, double *x)
{
  throw coupler::NotImplemented();
}


void BmiCoupler::
GetGridY(const int grid, double *y)
{
  throw coupler::NotImplemented();
}


void BmiCoupler::
GetGridZ(const int grid, double *z)
{
  throw coupler::NotImplemented();
}


int BmiCoupler::
GetGridNodeCount(const int grid)
{
  throw coupler::NotImplemented();
  /*
  if (grid == 0)
    return this->_model.shape[0];
  else
    return -1;
  */
}


int BmiCoupler::
GetGridEdgeCount(const int grid)
{
  throw coupler::NotImplemented();
}


int BmiCoupler::
GetGridFaceCount(const int grid)
{
  throw coupler::NotImplemented();
}


void BmiCoupler::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  throw coupler::NotImplemented();
}


void BmiCoupler::
GetGridFaceEdges(const int grid, int *face_edges)
{
  throw coupler::NotImplemented();
}


void BmiCoupler::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  throw coupler::NotImplemented();
}


void BmiCoupler::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  throw coupler::NotImplemented();
}


void BmiCoupler::
GetValue (std::string name, void *dest)
{
  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  nbytes = this->GetVarNbytes(name);
  
  memcpy (dest, src, nbytes);
}


void *BmiCoupler::
GetValuePtr (std::string name)
{
  if (name.compare("soil__storage") == 0)
    return (void*)this->_model.storage_m;
  else if (name.compare("soil__storage_change") == 0)
    return (void*)this->_model.storage_change_m;
  else  if (name.compare("soil__water_table") == 0)
    return (void*)this->_model.water_table_m;
  else if (name.compare("soil__moisture_content_total") == 0)
    return (void*)this->_model.SMCT;
  else if (name.compare("soil__moisture_content_layered") == 0)
    return (void*)this->_model.SMCL;
  else if (name.compare("soil__smc_profile_option_bmi") == 0)
    return (void*)this->_model.smc_profile_option_bmi;
  else {
    std::stringstream errMsg;
    errMsg << "variable "<< name << " does not exist";
    throw std::runtime_error(errMsg.str());
    return NULL;
  }
}


void BmiCoupler::
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


void BmiCoupler::
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


void BmiCoupler::
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


std::string BmiCoupler::
GetComponentName()
{
  return "Coupler BMI";
}


int BmiCoupler::
GetInputItemCount()
{
  return this->input_var_name_count;
}


int BmiCoupler::
GetOutputItemCount()
{
  return this->output_var_name_count;
}


std::vector<std::string> BmiCoupler::
GetInputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);

  return names;
}


std::vector<std::string> BmiCoupler::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->output_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


double BmiCoupler::
GetStartTime () {
  throw coupler::NotImplemented();
}


double BmiCoupler::
GetEndTime () {
  throw coupler::NotImplemented();
}


double BmiCoupler::
GetCurrentTime () {
  throw coupler::NotImplemented();
}


std::string BmiCoupler::
GetTimeUnits() {
  throw coupler::NotImplemented();
}


double BmiCoupler::
GetTimeStep () {
  throw coupler::NotImplemented();
}

#endif
