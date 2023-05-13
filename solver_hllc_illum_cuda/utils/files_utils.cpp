#if defined UTILS 

#include <iostream>
#include <string>

#if _MSC_VER
#define VERSION _MSVC_LANG
#else
#define VERSION __cplusplus
#endif

#if VERSION < 201703

#error "Need c++17 to use file_utils"

#else

#include <regex> //регул€рные выражени€
#include <filesystem> // need -std::c++17

namespace fs = std::filesystem;

static int delete_files(const std::string& path, const std::string& mask)
{
	const fs::path path_to_dir(path);
	const std::regex file_regex(mask);
	int count_del = 0;

	for (const auto& entry : fs::directory_iterator(path_to_dir)) //по всем файлам в директории
	{
		if (entry.is_regular_file())
		{
			std::string file_name = entry.path().filename().string();

			if (std::regex_match(file_name, file_regex))
			{
				printf("delete file: %s\n", entry.path().string().c_str());
				fs::remove(entry.path());
				count_del++;
			}
		}
	}

	return count_del;
}
static bool rename_files(const std::string& path, const std::string& mask, const int cur_idx, const int new_idx)
{
	const fs::path path_to_dir(path);
	const std::regex file_regex(mask);
	const std::regex replace_regex(std::to_string(cur_idx));

	bool exist_f = false;

	for (const auto& entry : fs::directory_iterator(path_to_dir)) //по всем файлам в директории
	{
		if (entry.is_regular_file())
		{
			std::string file_name = entry.path().filename().string();

			if (std::regex_match(file_name, file_regex))
			{
				std::string path = entry.path().generic_string();
				path.erase(path.end() - file_name.length(), path.end());// путь к файлу без имени

				std::string new_file_name = std::regex_replace(file_name, replace_regex, std::to_string(new_idx));
				fs::rename(entry.path(), path + new_file_name);

				exist_f = true;

				printf("rename: %s -> %s\n", file_name.c_str(), new_file_name.c_str());
			}
		}
	}

	return exist_f;
}
static int copy_files(const std::string& src_path, const std::string& dest_path, const std::string& mask)
{	
	if (!fs::exists(dest_path))
	{
		printf("no exist dest_dir: %s\n", dest_path.c_str());
		return 1;

		//fs::create_directory(dest_dir);
	}

	const fs::path path_to_dir(src_path);
	const std::regex file_regex(mask);
	int count_replace = 0;

	for (const auto& entry : fs::directory_iterator(path_to_dir)) //по всем файлам в директории
	{
		if (entry.is_regular_file())
		{
			std::string file_name = entry.path().filename().string();

			if (std::regex_match(file_name, file_regex))
			{
				fs::copy(entry.path(), dest_path, std::filesystem::copy_options::overwrite_existing);

				printf("copy file: %s\n", entry.path().string().c_str());
				count_replace++;
			}
		}
	}

	return count_replace;
}

int DeleteSolveFiles(int argc, char** argv)
{
	if (argc != 5)
	{
		printf("Error input data!\n");
		printf("Input:\n");
		printf("path\n");
		printf("main_name_file (like \"Solve\")\n");
		printf("max number of files idx\n");
		printf("del mask\n");
		return 1;
	}

	const std::string path = argv[1];
	const std::string main_name = argv[2];
	const int n_files = std::stoi(argv[3]);
	const int del_mask = std::stoi(argv[4]);

	for (int i = 0; i < n_files; i++) //первый сохран€ем
	{
		if (i % del_mask != 0) // оставить только каждый del_mask file (ex: 0,2,4,6,...)
		{
			std::string mask_file_regex(main_name + std::to_string(i) + "\[^0-9]*"); //все кроме цифр
			delete_files(path, mask_file_regex);
		}
	}

	return 0;
}
int ReduceNameSolveFiles(int argc, char** argv)
{
	if (argc != 4)
	{
		printf("Error input data!\n");
		printf("Input:\n");
		printf("path\n");
		printf("main_name_file (like \"Solve\")\n");
		printf("max number of files idx\n");
		return 1;
	}

	const std::string path = argv[1];
	const std::string main_name = argv[2];
	const int n_files = std::stoi(argv[3]);		

	int new_idx = 0;
	for (int i = 0; i < n_files; i++)
	{
		std::string mask_file_regex(main_name + std::to_string(i) + "\[^0-9]*");

		if (rename_files(path, mask_file_regex, i, new_idx))
		{
			new_idx++;
		}
	}
}
int CopySolveFiles(int argc, char** argv)
{
	if (argc != 6)
	{
		printf("Error input data!\n");
		printf("Input:\n");
		printf("src_path\n");
		printf("dest_path\n");
		printf("main_name_file (like \"Solve\")\n");
		printf("max number of files idx\n");
		printf("del mask\n");
		return 1;
	}

	const std::string src_path = argv[1];
	const std::string dest_path = argv[2];
	const std::string main_name = argv[3];
	const int n_files = std::stoi(argv[4]);
	const int cpy_mask = std::stoi(argv[5]);

	for (int i = 0; i < n_files; i += cpy_mask)
	{
		std::string mask_file_regex(main_name + std::to_string(i) + "\[^0-9]*");
		copy_files(src_path, dest_path, mask_file_regex);
	}
}

#include <fstream>
#include <iterator>
#include "../global_def.h"
int MakeSimpleTxtFromData(int argc, char** argv)
{
	if (argc != 3)
	{
		printf("Error input data!\n");
		printf("Input:\n");
		printf("path\\data.txt\n");
		printf("path\\out_simple_format.txt\n");
		return 1;
	}

	std::string file_data = argv[1];
	std::string file_simple_format = argv[2];

	//std::string file_data = "D:\\Desktop\\FilesCourse\\Test\\auto_metis\\cone_m_epart_2.txt";
	//std::string file_simple_format = "D:\\Desktop\\FilesCourse\\Test\\auto_metis\\cone_m_epart_2_simple.txt";


	std::ifstream ifile;
	OPEN_FSTREAM(ifile, file_data.c_str());

	std::string file_content; // ((std::istream_iterator<char>(ifile.rdbuf())), std::istream_iterator<char>());
	file_content.assign((std::istreambuf_iterator<char>(ifile.rdbuf())), std::istreambuf_iterator<char>());
	ifile.close();

	OPEN_FSTREAM(ifile, file_data.c_str());
	int n = std::count(std::istreambuf_iterator<char>(ifile), std::istreambuf_iterator<char>(), '\n');
	ifile.close();

	file_content = std::to_string(n) + "\n" + file_content;

	std::ofstream ofile;
	OPEN_FSTREAM(ofile, file_simple_format.c_str());
	ofile << file_content;
	ofile.close();

	return 0;
}
#endif
#endif //UTILS
