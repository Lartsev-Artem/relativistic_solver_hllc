#include "illum_utils.h"
#if defined ILLUM

#include "../../file_module/reader_bin.h"

/*! Пересчет излучения с грани в ячейку в заданном направлении
num_dir - номер направления
illum_on_face - излучение на гранях по всем направлениям
*/
int GetDirectionIllumFromFace(const int size_grid, const int num_dir, const std::vector<Type>& illum_on_face, std::vector<Type>& illum_in_cell)
{	
	if (illum_on_face.size() < size_grid* base +size_grid * num_dir * base) RETURN_ERR("illum_on_face hasn't enough data\n");

	illum_in_cell.resize(size_grid, 0);	
	for (size_t j = 0; j < size_grid; j++)
	{
		const int N = num_dir * base * size_grid + j * base;
		for (size_t k = 0; k < base; k++)
		{
			illum_in_cell[j] += illum_on_face[N + k];
		}
		illum_in_cell[j] /=  base;		
	}
	return 0;
}
#endif
#if defined ILLUM && defined SOLVE
int SolveIllumAndHLLC(const Type tau, std::vector<elem_t>& cells)
{
	for (auto& el : cells)
	{
		//el.phys_val.d += 0;
		//el.phys_val.v += (C_LIGHT*DIST/TIME/PRESSURE) * tau * (-el.illum_val.div_impuls * C_LIGHT_INV);  //- tau*(stream[i][0] - prev_stream[i][0]) / tau / c / c);  // +F_g //vx
		//el.phys_val.p += (1./TIME/RADIATION)*tau * (-el.illum_val.div_stream); //  tau*(-(energy[i] - prev_energy[i]) / tau / c)  //+F_g.dot(vel)  //e

		el.conv_val.v +=   tau * (-el.illum_val.div_impuls * 1);  //- tau*(stream[i][0] - prev_stream[i][0]) / tau / c / c);  // +F_g //vx
		el.conv_val.p +=   tau * (-el.illum_val.div_stream); //  tau*(-(energy[i] - prev_energy[i]) / tau / c)  //+F_g.dot(vel)  //e
	}	
	return 0;
}

int InitIllum(file_name main_dir, grid_t& grid)
{
	double _clock = -omp_get_wtime();

	if (ReadDataArray(solve_mode, main_dir, grid))
	{
		RETURN_ERR("Error reading the data array vtk\n");
	}
	_clock += omp_get_wtime();

	WRITE_LOG("Reading the data array " << _clock << " c\n");	

#if 0 // инициализация составляющий для излучения (не надо для квазилинейного случая)
	CalculateIllumOptMemory(main_dir, use_cuda, cuda_mod, accuracy, max_number_of_iter, directions, Illum, int_scattering, squares,
		grid, neighbours_id_faces, OutC, Out, In, S, X, X0,
		res_inner_bound, id_try_surface, sorted_id_cell, ShiftOut, ShiftRes, ShiftX0, ShiftTry, U_full_prev);

	CalculateIllumParamOptMemory(size_grid, count_directions,
		directions, squares, normals, squares_cell, volume,
		Illum, energy, stream, impuls, div_stream, div_impuls);

	printf("Init Illum completed\n\n");
#endif

	return 0;
}

#endif //ILLUM


Matrix3 IntegarteDirection9(const int num_cell, const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface) {
	Matrix3 res = Matrix3::Zero();
	int n = squares.size();  // number of directions
	int m = Illum.size() / n;  // number of cells

	for (size_t i = 0; i < 3; i++)
		for (size_t k = 0; k < 3; k++)

			for (size_t j = 0; j < n; j++) {
				res(i, k) += directions[j][i] * directions[j][k] * (Illum[m * j + num_cell] * squares[j]);
			}

	return res / scuare_surface;
}
int MakeImpuls(const vector<Type>& Illum, const vector<Vector3>& directions, const vector<Type>& squares, const Type scuare_surface, vector<Matrix3>& impuls) {

	const int n = impuls.size();

	for (size_t i = 0; i < n; ++i) {
		impuls[i] = IntegarteDirection9(i, Illum, directions, squares, scuare_surface);
	}

	return 0;
}
