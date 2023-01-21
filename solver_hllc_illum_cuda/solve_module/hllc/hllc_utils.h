#ifndef HLLC_UTILS_H
#define HLLC_UTILS_H

#include "../solve_config.h"

#ifdef HLLC

#include "../solve_global_struct.h"

#include "../../global_def.h"
#include "../../global_value.h"

int HllcPhysToConv(std::vector<elem_t>& cells);
int HllcConvToPhys(std::vector<elem_t>& cells);
int HllcGetTimeStep(hllc_value_t& hllc_set, const std::vector<elem_t>& cells);

int InitHLLC(file_name file_settings_hllc, hllc_value_t& hllc_set, file_name file_init_value, std::vector<elem_t>& cells);

int HLLC_3d(const Type tau, grid_t& grid);

#endif // HLLC
#endif //HLLC_UTILS_H