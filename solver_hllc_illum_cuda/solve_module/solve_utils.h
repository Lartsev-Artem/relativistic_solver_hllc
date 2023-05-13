#if !defined SOLVE_UTILS_H && defined SOLVE
#define SOLVE_UTILS_H
#include "solve_config.h"

#include "solve_global_struct.h"

int HLLC_STEP(const Type tau, grid_t& grid);
int GetTimeStep(hllc_value_t& hllc_cfg, const grid_t& grid);
int HLLC_INIT(file_name file_settings_hllc, hllc_value_t& hllc_set, file_name file_init_value, std::vector<elem_t>& cells);

void Init_MPI();

int StartLowDimensionTask(file_name main_dir);

int TestDivStream(const std::vector<Vector3>& centers_face, grid_t& grid);
int TestDivStream(file_name base_adress);

#ifdef USE_MPI

int GetSend(const int myid);
int GetDisp(const int myid);

void SendPhysValue(flux_t* phys, const int size, const int msg_id);
int MPI_RHLLC_3d(const int myid, const Type tau, grid_t& grid);
int RHLLC_3d_MPI(const Type tau, grid_t& grid);

void GetSend(const int np, const int n, std::vector<int>& send_count);
void GetDisp(const int np, const int n, std::vector<int>& disp);
void GetDispSend(const int np, const int n, const int coef, std::vector<int>& send_count, std::vector<int>& disp);
#endif //USE_MPI

#endif //SOLVE_UTILS_H