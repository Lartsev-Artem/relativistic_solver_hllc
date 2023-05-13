#ifndef RHLLC_UTILS_H
#define RHLLC_UTILS_H

#include "../solve_config.h"

#ifdef RHLLC

#include "../solve_global_struct.h"

#include "../../global_def.h"
#include "../../global_value.h"

int RHllcPhysToConv(std::vector<elem_t>& cells);

int ReBuildPhysicValue(const Vector4& U, Vector4& W);
int ReBuildPhysicValue(const  std::vector<Vector4>& U, std::vector<Vector4>& W);
int ReBuildConvValue(const std::vector<Vector4>& W, std::vector<Vector4>& U);

int rhllc_get_conv_value_ost1098(const flux_t& W, flux_t& U);

int InitRHLLC(file_name file_settings_hllc, hllc_value_t& hllc_set, file_name file_init_value, std::vector<elem_t>& cells);

int RHLLC_3d(const Type tau, grid_t& grid);

int RHllcGetTimeStep(hllc_value_t& hllc_set, const std::vector<elem_t>& cells);

Type Density(const Vector3& p);
Type Pressure(const Vector3& p);
Vector3 Velocity(const Vector3& p);

#if NUMBER_OF_MEASUREMENTS == 3 && defined RHLLC_MPI
int ReadMpiConf(const std::string& file, int myid, int np);
int InitMPI_RHllc(const std::vector<elem_t>& cells);
#endif

#endif //RHLLC
#endif //RHLLC_UTILS_H