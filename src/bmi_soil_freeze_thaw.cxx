#ifndef BMI_SFT_C_INCLUDED
#define BMI_SFT_C_INCLUDED

#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "../include/bmi_soil_freeze_thaw.hxx"
#include "../include/soil_freeze_thaw.hxx"
#include "../include/Logger.hpp"
#include "../bmi/bmi.hxx"
#include <algorithm>

#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

BmiSoilFreezeThaw::BmiSoilFreezeThaw() : m_serialized_vec{} {
  this->input_var_names[0]  = "ground_temperature";
  this->input_var_names[1]  = "soil_moisture_profile";
 
  this->output_var_names[0] = "ice_fraction_schaake";
  this->output_var_names[1] = "ice_fraction_xinanjiang";
  this->output_var_names[2] = "num_cells";
  this->output_var_names[3] = "soil_temperature_profile";
  this->output_var_names[4] = "soil_ice_fraction";
  this->output_var_names[5] = "ground_heat_flux";
 
  // add calibratable parameters
  this->calib_var_names[0]  = "smcmax";
  this->calib_var_names[1]  = "b";
  this->calib_var_names[2]  = "satpsi";
 
  // ensure empty serialized state
  this->m_serialized_length = 0;
};

void BmiSoilFreezeThaw::
Initialize (std::string config_file)
{
  // Initialize the Error, Warning and Trapping System
#ifdef EWTS_HAVE_NGEN_BRIDGE    
  EwtsInit(SFT_MODULE_ID, true);
#else
  EwtsInit(SFT_MODULE_ID, false);
#endif
  LOG(LogLevel::INFO, "Initializing SFT");

  if (config_file.compare("") != 0 )
      this->state = new soilfreezethaw::SoilFreezeThaw(config_file);

    verbosity= this->state->verbosity;
}

void BmiSoilFreezeThaw::
Update()
{
  this->state->Advance();
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
    this->state->dt = frac * dt;
    this->state->Advance();
    this->state->dt = dt;
  }
}


void BmiSoilFreezeThaw::
Finalize()
{
  if (this->state != nullptr)
    delete this->state;
}

int BmiSoilFreezeThaw::
GetVarGrid(std::string name)
{
  if (name.compare("num_cells") == 0 
     || name.compare("ice_fraction_scheme_bmi") == 0
     || name.compare("serialization_free") == 0)
    return 0; // int
  else if (name.compare("ground_temperature") == 0 || name.compare("ice_fraction_schaake") == 0
	   || name.compare("ice_fraction_xinanjiang") == 0 || name.compare("soil_ice_fraction") == 0
	   || name.compare("ground_heat_flux") == 0 || name.compare("smcmax") == 0
	   || name.compare("b") == 0 || name.compare("satpsi") == 0
     || name == "reset_time")
    return 1; //double
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_temperature_profile") == 0)
    return 2; // arrays
  else if (name.compare("serialization_state") == 0)
    return 3; // char
  else if (name.compare("serialization_create") == 0
           || name.compare("serialization_size") == 0)
    return 4; // unit64_t
  else
    return -1;
}


std::string BmiSoilFreezeThaw::
GetVarType(std::string name)
{
  int grid_id = GetVarGrid(name);

  if (grid_id == 0)
    return "int";
  else if (grid_id == 1 || grid_id == 2)
    return "double";
  else if (grid_id == 3)
    return "char";
  else if (grid_id == 4)
    return "uint64_t";
  else
    return "";
}


int BmiSoilFreezeThaw::
GetVarItemsize(std::string name)
{
  int grid_id = GetVarGrid(name);

  if (grid_id == 0)
    return sizeof(int);
  else if (grid_id == 1 || grid_id == 2)
    return sizeof(double);
  else if (grid_id == 3)
    return sizeof(char);
  else if (grid_id == 4)
    return sizeof(uint64_t);
  else
    return 0;
  
}


std::string BmiSoilFreezeThaw::
GetVarUnits(std::string name)
{
  if (name.compare("ground_temperature") == 0 || name.compare("soil_temperature_profile") == 0)
    return "K";
  else if (name.compare("ground_heat_flux") == 0)
    return "W m-2";
  else if (name.compare("ice_fraction_schaake") == 0 ||
           name.compare("ice_fraction_xinanjiang") == 0 ||
           name.compare("soil_ice_fraction") == 0 ||
           name.compare("soil_moisture_profile") == 0)
    return "1"; // UDUNITS dimensionless
  else
    return "none";
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
  else if (name.compare("ice_fraction_xinanjiang") == 0 || name.compare("soil_ice_fraction") == 0)
    return "node";
  else if (name.compare("ice_fraction_schaake") == 0 ||  name.compare("num_cells") == 0
	   || name.compare("ground_heat_flux") == 0)
    return "node";
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_temperature_profile") == 0)
    return "node";
  else
    return "";
}


void BmiSoilFreezeThaw::
GetGridShape(const int grid, int *shape)
{
  if (grid == 2) {
    shape[0] = this->state->shape[0];
  }
}


void BmiSoilFreezeThaw::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->state->spacing[0];
  }
}


void BmiSoilFreezeThaw::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->state->origin[0];
  }
}


int BmiSoilFreezeThaw::
GetGridRank(const int grid)
{
  if (grid == 0 || grid == 1 || grid == 2)
    return 1;
  else
    return -1;
}


int BmiSoilFreezeThaw::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 1) // for scalars
    return 1;
  else if (grid == 2)        // for arrays
    return this->state->shape[0];
  else if (grid == 3)        // serialized data
    return this->m_serialized_length;
  else if (grid == 4)
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
  LOG(LogLevel::WARNING,"Not implemented in SFT");
  throw std::logic_error("Not Implemented in SFT");
}


void BmiSoilFreezeThaw::
GetGridY(const int grid, double *y)
{
  LOG(LogLevel::WARNING,"Not implemented in SFT");
  throw std::logic_error("Not Implemented in SFT");
}


void BmiSoilFreezeThaw::
GetGridZ(const int grid, double *z)
{
  LOG(LogLevel::WARNING,"Not implemented in SFT");
  throw std::logic_error("Not Implemented in SFT");
}


int BmiSoilFreezeThaw::
GetGridNodeCount(const int grid)
{
  if (grid == 0)
    return this->state->shape[0];
  else
    return -1;
}


void BmiSoilFreezeThaw::
GetValue (std::string name, void *dest)
{
  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  
  if (name.compare("serialization_state") == 0) {
    memcpy(dest, src, this->m_serialized_length);
  } else {
    nbytes = this->GetVarNbytes(name);
    memcpy (dest, src, nbytes);
  }
}


void *BmiSoilFreezeThaw::
GetValuePtr (std::string name)
{
  if (name.compare("soil_temperature_profile") == 0)
    return (void*)this->state->soil_temperature.data();
  if (name.compare("soil_moisture_profile") == 0)
    return (void*)this->state->soil_moisture_content.data();
  else if (name.compare("ground_temperature") == 0 )
    return (void*)(&this->state->ground_temp);
  else if (name.compare("ground_heat_flux") == 0)
    return (void*)(&this->state->ground_heat_flux);
  else if (name.compare("num_cells") == 0)
    return (void*)(&this->state->ncells);
  else if (name.compare("ice_fraction_schaake") == 0)
    return (void*)(&this->state->ice_fraction_schaake);
  else if (name.compare("ice_fraction_xinanjiang") == 0)
    return (void*)(&this->state->ice_fraction_xinanjiang);
  else if (name.compare("soil_ice_fraction") == 0)
    return (void*)(&this->state->soil_ice_fraction);
  else if (name.compare("ice_fraction_scheme_bmi") == 0)
    return (void*)(&this->state->ice_fraction_scheme_bmi); // leaving this for now, but probably not needed/used
  else if (name.compare("smcmax") == 0)
    return (void*)(&this->state->smcmax);
  else if (name.compare("b") == 0)
    return (void*)(&this->state->b);
  else if (name.compare("satpsi") == 0)
    return (void*)(&this->state->satpsi);
  else if (name.compare("serialization_state") == 0)
    return (void*)(this->m_serialized_vec.data());
  else if (name.compare("serialization_size") == 0) {
    return (void*)(&this->m_serialized_length);
  }
  else {
    //std::stringstream errMsg;
    //errMsg << "variable "<< name << " does not exist";
    std::string errMsg = "Variable " + name + " does not exist\n";
    LOG(LogLevel::WARNING, errMsg);
    throw std::runtime_error(errMsg);
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

  // special case for clearing serialized data
  if (name.compare("serialization_free") == 0) {
    this->clear_serialized();
    return;
  } else if (name.compare("serialization_state") == 0) {
    this->load_serialized((char*)src);
    return;
  } else if (name.compare("serialization_create") == 0) {
    this->new_serialized();
    return;
  } else if (name == "reset_time") {
    this->reset_time();
    return;
  } else {
    dest = this->GetValuePtr(name);
  }
  
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
  
  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);
  
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
  return this->state->endtime;
}


double BmiSoilFreezeThaw::
GetCurrentTime () {
  return this->state->time;
}


std::string BmiSoilFreezeThaw::
GetTimeUnits() {
  return "s";
}


double BmiSoilFreezeThaw::
GetTimeStep () {
  return this->state->dt;
}


int BmiSoilFreezeThaw::
GetGridEdgeCount(const int grid)
{
  LOG(LogLevel::WARNING,"Not implemented in SFT");
  throw std::logic_error("Not Implemented in SFT");

  return 0;
}


int BmiSoilFreezeThaw::
GetGridFaceCount(const int grid)
{
  LOG(LogLevel::WARNING,"Not implemented in SFT");
  throw std::logic_error("Not Implemented in SFT");

  return 0;
}


void BmiSoilFreezeThaw::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  LOG(LogLevel::WARNING,"Not implemented in SFT");
  throw std::logic_error("Not Implemented in SFT");
}


void BmiSoilFreezeThaw::
GetGridFaceEdges(const int grid, int *face_edges)
{
  LOG(LogLevel::WARNING,"Not implemented in SFT");
  throw std::logic_error("Not Implemented in SFT");
}


void BmiSoilFreezeThaw::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  LOG(LogLevel::WARNING,"Not implemented in SFT");
  throw std::logic_error("Not Implemented in SFT");
}


void BmiSoilFreezeThaw::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  LOG(LogLevel::WARNING,"Not implemented in SFT");
  throw std::logic_error("Not Implemented in SFT");
}


template<class Archive>
void BmiSoilFreezeThaw::
serialize(Archive& ar, const unsigned int version) {
  // could throw archive_exception if the size of the archive array isn't the same as the size of the BMI
  soilfreezethaw::SoilFreezeThaw* state = this->state;
  // in Advance
  ar & state->time;
  ar & state->soil_temperature_prev;
  ar & state->soil_liquid_content;

  // in ThermalConductivity
  ar & state->thermal_conductivity;

  // in SoilHeatCapacity
  ar & state->heat_capacity;

  // in SolveDiffusionsEquation
  ar & state->ground_heat_flux;
  ar & state->bottom_heat_flux;
  ar & state->soil_temperature;

  // in PhaseChange
  // soil_temperature; already covered above
  ar & state->energy_consumed;
  // ar & make_array(state->soil_liquid_content, size);
  ar & state->soil_moisture_content;
  ar & state->soil_ice_content;

  // in ComputeIceFraction
  // soil_ice_content; already covered above
  ar & state->ice_fraction_schaake;
  ar & state->ice_fraction_xinanjiang;
  ar & state->soil_ice_fraction;
  ar & state->ice_fraction_scheme_bmi;

  // in EnergyBalanceCheck
  ar & state->energy_balance;
}


void BmiSoilFreezeThaw::
new_serialized() {
  // resize to reserve space for the serialized size as a header
  this->m_serialized_vec.resize(sizeof(uint64_t));
  boost::archive::binary_oarchive archive(this->m_serialized_vec);
  try {
    archive << (*this);
    this->m_serialized_length = this->m_serialized_vec.size();
    // copy size of serialized data to the beginning as a header
    uint64_t serialized_size = this->m_serialized_length - sizeof(uint64_t);
    memcpy(this->m_serialized_vec.data(), &serialized_size, sizeof(uint64_t));
  } catch (const std::exception &e) {
    LOG(LogLevel::SEVERE, "Serializing SFT encountered an error: %s", e.what());
    this->m_serialized_length = 0;
    throw;
  }
}


void BmiSoilFreezeThaw::
load_serialized(char* data) {
  // pull data size from header of raw data ptr
  uint64_t size;
  memcpy(&size, data, sizeof(uint64_t));
  // create stream from after the size header
  membuf stream(data + sizeof(uint64_t), size);
  boost::archive::binary_iarchive archive(stream);
  try {
    archive >> (*this);
    this->clear_serialized();
  } catch (const std::exception &e) {
    LOG(LogLevel::SEVERE, "Deserializing SFT encounterd an error: %s", e.what());
    throw;
  }
}


// Clear the currently saved serialized data and release memory
void BmiSoilFreezeThaw::
clear_serialized() {
  this->m_serialized_vec.clear();
  this->m_serialized_vec.shrink_to_fit();
  this->m_serialized_length = 0;
}

void BmiSoilFreezeThaw::
reset_time() {
  // time doesn't seem to be used anywhere but GetCurrentTime, so safe to set to 0 and nothing else
  this->state->time = 0.0;
}


#endif
