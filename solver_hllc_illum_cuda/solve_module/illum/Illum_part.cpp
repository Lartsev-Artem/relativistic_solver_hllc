#include "solve_short_characteristics_main.h"

#ifdef USE_VTK
#else

int CalculateIllum(const std::string& main_dir, bool& use_cuda, const int cuda_mod, const Type accuracy, const int max_number_of_iter,
	const std::vector<Vector3>& directions, std::vector<Type>& Illum, std::vector<Type>& int_scattering, std::vector<Type> squares,
	std::vector<cell>& grid, const std::vector<int>& neighbours_id_faces, const std::vector<ShortId>& OutC,
	const std::vector<ShortId>& Out, const std::vector<ShortId>& In, const std::vector<Type>& S, const std::vector<Vector3>& X, const std::vector<Vector2>& X0,
	const std::vector<Type>& res_inner_bound, const std::vector<int>& id_try_surface, const vector<int>& sorted_id_cell,
	const vector<uint64_t>& ShiftOut, const vector<uint64_t>& ShiftRes, const vector<uint64_t>& ShiftX0, const vector<int>& ShiftTry,
	const std::vector<VectorX>& U_full)
{
	
	const int count_directions = directions.size();
	const int count_cells = size_grid;

	// вроде не об€зательно. ѕ–ќ¬≈–»“№
	//Illum.assign(4 * count_cells * count_directions, 0);
	//int_scattering.assign(count_cells * count_directions, 0);

	int count = 0;
	Type norm = 0;

#ifdef WRITE_LOG
	ofstream ofile;
	ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
#endif

	do {
		Type _clock = -omp_get_wtime();

		norm = -1;
		/*---------------------------------- далее FOR по направлени€м----------------------------------*/

		for (register int num_direction = 0; num_direction < count_directions; ++num_direction)
		{
			int posX = 0;
			int posX0 = 0;
			int posOut = 0;
			int posIn = 0;
			int posS = 0;

			int pos_in_res = 0;
			int id_try_pos = 0;


			register int num_cell;
			ShortId n_out;
			Vector3 I;
			Type sumI = 0;


			/*---------------------------------- далее FOR по €чейкам----------------------------------*/
			for (int h = 0; h < count_cells; ++h) {
				num_cell = sorted_id_cell[num_direction * count_cells + h];

				n_out = OutC[num_direction * count_cells + h];

				sumI = 0;

				for (ShortId i = 0; i < n_out; ++i) 
				{
					ShortId num_out_face = Out[ShiftOut[num_direction] + posOut++];

					//GetNodes
					for (int num_node = 0; num_node < 3; ++num_node) {
						Vector3 x = X[3 * ShiftOut[num_direction] + posX++];

						ShortId num_in_face = In[3 * ShiftOut[num_direction] + posIn++];

						Type I_x0 = CalculateIllumeOnInnerFace(num_cell, num_direction, num_in_face, x, X0, grid, neighbours_id_faces,
							id_try_pos, pos_in_res, posX0, ShiftRes[num_direction], ShiftX0[num_direction], ShiftTry[num_direction]);

						Type s = S[3 * ShiftOut[num_direction] + posS++];

						I[num_node] = CurGetIllum(num_cell, num_direction, x, s, I_x0, int_scattering, U_full);
						sumI += I[num_node];

					}//num_node

					// если хранить не знгачени€ у коэфф. интерпол€ции
					{
						//Vector3 coef;
						//if (num_out_face == 3)
						//	coef = inclined_face_inverse * I;// GetInterpolationCoefInverse(inclined_face_inverse, I);
						//else
						//	coef = straight_face_inverse * I;// GetInterpolationCoefInverse(straight_face_inverse, I);
					}

					grid[num_cell].nodes_value[num_out_face] = I; //coef;  // храним значение узлов а не коэф. интерпол€ции

					//Illum[num_direction * (4*count_cells) + num_cell * 4 + num_out_face] = (I[0] + I[1] + I[2]) / 3;					

					int neighbor_id_face = neighbours_id_faces[num_cell * 4 + num_out_face];
					if (neighbor_id_face >= 0)
						grid[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] = I;// coef;


				} //num_out_face

				for (size_t i = 0; i < 4; i++)
				{
					Vector3 Il = grid[num_cell].nodes_value[i];
					const Type curI = (Il[0] + Il[1] + Il[2]) / 3;
					const int id = num_direction * (4 * count_cells) + num_cell * 4 + i;
					{
						Type buf_norm = fabs(Illum[id] - curI) / curI;
						if (buf_norm > norm) norm = buf_norm;
					}
					Illum[id] = curI;
				}
				//illum[num_direction * count_cells + num_cell] = sumI / (3 * n_out);				
			}
			/*---------------------------------- конец FOR по €чейкам----------------------------------*/

			//std::cout << "End direction number: " << num_direction << '\n';
		}
		/*---------------------------------- конец FOR по направлени€м----------------------------------*/


	
		if (max_number_of_iter > 1)  // пропуск первой итерации
		{
			if (use_cuda)
			{
				if (CalculateIntScattering(32, count_cells, count_directions, Illum, int_scattering)) { // если была ошибка пересчитать на CPU 
					CalculateIntOmp(count_cells, count_directions, Illum, directions, squares, int_scattering);
					ClearDevice(cuda_mod);
					use_cuda = false;
				}
			}
			else
				CalculateIntOmp(count_cells, count_directions, Illum, directions, squares, int_scattering);
		}

		//Illum.swap(Illum2);

		_clock += omp_get_wtime();
		
#ifdef WRITE_LOG
		printf("Error: %lf\n", norm);
		printf("Time of iter: %lf\n", _clock);
		printf("End iter_count number: %d\n", count);
		
		ofile << "Error:= " << norm << '\n';
		ofile << "Time of iter: " << _clock << '\n';
		ofile << "End iter_count number: " << count << '\n';
#endif
		count++;
		
	} while (norm > accuracy && count < max_number_of_iter);

#ifdef WRITE_LOG
	ofile.close();
#endif

	return 0;
}


int CalculateIllumOptMemory(const std::string& main_dir, bool& use_cuda, const int cuda_mod, const Type accuracy, const int max_number_of_iter,
	const std::vector<Vector3>& directions, std::vector<Type>& Illum, std::vector<Type>& int_scattering, std::vector<Type> squares,
	
	std::vector<cell>& grid, const std::vector<int>& neighbours_id_faces, const std::vector<ShortId>& OutC,
	const std::vector<ShortId>& Out, const std::vector<ShortId>& In,  std::vector<Type>& S,  std::vector<Vector3>& X,  std::vector<Vector2>& X0,
	const std::vector<Type>& res_inner_bound, const std::vector<int>& id_try_surface, const vector<int>& sorted_id_cell,
	const vector<uint64_t>& ShiftOut, const vector<uint64_t>& ShiftRes, const vector<uint64_t>& ShiftX0, const vector<int>& ShiftTry,
	const std::vector<VectorX>& U_full)
{

	std::string name_file_local_x0 = main_dir + "LocX0";
	std::string name_file_x = main_dir + "X";
	std::string name_file_s = main_dir + "S";
	std::string name_file_sizes = main_dir + "Size";	
	
	const int count_directions = directions.size();
	const int count_cells = size_grid;

	int count = 0;
	Type norm = 0;
	Illum.resize(4*size_grid, 0);

#ifdef WRITE_LOG
	ofstream ofile;
	//ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
#endif
	int countX, countX0, countOutC, countOut, countIn, countS, countRes, countTry;
	ReadSizes(name_file_sizes, countX, countX0, countOutC, countOut, countIn, countS, countRes, countTry);
	X0.clear();
	do {
		Type _clock = -omp_get_wtime();

		norm = -1;
		/*---------------------------------- далее FOR по направлени€м----------------------------------*/
		FILE* fs;
		FILE* fx;
		FILE* fI;
		fs = fopen((name_file_s + ".bin").c_str(), "rb");
		fx = fopen((name_file_x + ".bin").c_str(), "rb");
		fI = fopen((main_dir + "FormIllum.bin").c_str(), "wb");


		for (register int num_direction = 0; num_direction < count_directions; ++num_direction)
		{
#ifdef WRITE_LOG
			ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
			ofile << "Iter direction: " << num_direction << '\n';
			ofile.close();
#endif
			Illum.assign(4*size_grid, 0);
			if (num_direction != count_directions - 1) {
				int n = 3*ShiftOut[num_direction + 1];
				S.resize(n);
				X.resize(n);
			}
			else 
			{
				int n = countS - 3 * ShiftOut[num_direction]; 		//последнее направление
				S.resize(n);
				X.resize(n);		
			}
		
			fread_unlocked(S.data(), sizeof(Type), S.size(), fs);									
			fread(X.data(), sizeof(Vector3), X.size(), fx);
						

		/*	f = fopen((name_file_local_x0 + ".bin").c_str(), "rb");
			if (!f) { printf("Err x0\n"); return 1; }
			fread(X0.data(), sizeof(Vector2), X0.size(), f);
			fclose(f);
			printf("Read x0 files\n");*/


			int posX = 0;
			int posX0 = 0;
			int posOut = 0;
			int posIn = 0;
			int posS = 0;

			int pos_in_res = 0;
			int id_try_pos = 0;


			register int num_cell;
			ShortId n_out;
			Vector3 I;
			Type sumI = 0;


			/*---------------------------------- далее FOR по €чейкам----------------------------------*/
			for (int h = 0; h < count_cells; ++h) {
				num_cell = sorted_id_cell[num_direction * count_cells + h];

				n_out = OutC[num_direction * count_cells + h];

				sumI = 0;

				for (ShortId i = 0; i < n_out; ++i)
				{
					ShortId num_out_face = Out[ShiftOut[num_direction] + posOut++];

					//GetNodes
					for (int num_node = 0; num_node < 3; ++num_node) {
						Vector3 x = X[posX++];

						ShortId num_in_face = In[3 * ShiftOut[num_direction] + posIn++];

						Type I_x0 = CalculateIllumeOnInnerFace(num_cell, num_direction, num_in_face, x, X0, grid, neighbours_id_faces,
							id_try_pos, pos_in_res, posX0, ShiftRes[num_direction], ShiftX0[num_direction], ShiftTry[num_direction]);

						Type s = S[posS++];

						I[num_node] = CurGetIllum(num_cell, num_direction, x, s, I_x0, int_scattering, U_full);
						sumI += I[num_node];

					}//num_node

					// если хранить не знгачени€ у коэфф. интерпол€ции
					{
						//Vector3 coef;
						//if (num_out_face == 3)
						//	coef = inclined_face_inverse * I;// GetInterpolationCoefInverse(inclined_face_inverse, I);
						//else
						//	coef = straight_face_inverse * I;// GetInterpolationCoefInverse(straight_face_inverse, I);
					}

					grid[num_cell].nodes_value[num_out_face] = I; //coef;  // храним значение узлов а не коэф. интерпол€ции

					//Illum[num_direction * (4*count_cells) + num_cell * 4 + num_out_face] = (I[0] + I[1] + I[2]) / 3;					

					int neighbor_id_face = neighbours_id_faces[num_cell * 4 + num_out_face];
					if (neighbor_id_face >= 0)
						grid[neighbor_id_face / 4].nodes_value[neighbor_id_face % 4] = I;// coef;


				} //num_out_face

				for (size_t i = 0; i < 4; i++)
				{
					Vector3 Il = grid[num_cell].nodes_value[i];
					const Type curI = (Il[0] + Il[1] + Il[2]) / 3;
					//const int id = num_direction * (4 * count_cells) + num_cell * 4 + i;
					const int id = num_cell * 4 + i;
					{
						Type buf_norm = fabs(Illum[id] - curI) / curI;
						if (buf_norm > norm) norm = buf_norm;
					}
					Illum[id] = curI;
				}
				//illum[num_direction * count_cells + num_cell] = sumI / (3 * n_out);				
			}
			/*---------------------------------- конец FOR по €чейкам----------------------------------*/
			{								
					for (size_t j = 0; j < count_cells; j++)
					{
						const int N = j * 4;
						Type I = (Illum[N] + Illum[N + 1] + Illum[N + 2] + Illum[N + 3]) / 4;
						fwrite(&I, sizeof(Type), 1, fI);						
					}
			}
			
			//std::cout << "End direction number: " << num_direction << '\n';
		}
		/*---------------------------------- конец FOR по направлени€м----------------------------------*/

		fclose(fs);
		fclose(fx);
		fclose(fI);

		if (max_number_of_iter > 1)  // пропуск первой итерации
		{
			if (use_cuda)
			{
				if (CalculateIntScattering(32, count_cells, count_directions, Illum, int_scattering)) { // если была ошибка пересчитать на CPU 
					CalculateIntOmp(count_cells, count_directions, Illum, directions, squares, int_scattering);
					ClearDevice(cuda_mod);
					use_cuda = false;
				}
			}
			else
				CalculateIntOmp(count_cells, count_directions, Illum, directions, squares, int_scattering);
		}

		//Illum.swap(Illum2);

		_clock += omp_get_wtime();

#ifdef WRITE_LOG
		printf("Error: %lf\n", norm);
		printf("Time of iter: %lf\n", _clock);
		printf("End iter_count number: %d\n", count);

		ofile.open(main_dir + "File_with_Logs_solve.txt", std::ios::app);
		ofile << "Error:= " << norm << '\n';
		ofile << "Time of iter: " << _clock << '\n';
		ofile << "End iter_count number: " << count << '\n';
		ofile.close();
#endif
		count++;

	} while (norm > accuracy && count < max_number_of_iter);


	return 0;
}
#endif // USE_VTK