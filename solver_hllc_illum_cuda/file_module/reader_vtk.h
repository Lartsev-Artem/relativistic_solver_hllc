#include "../prj_config.h"

#if (!defined READER_VTK && defined USE_VTK)

#define READER_VTK

#include "../global_def.h"

#include<vtk-9.0/vtkUnstructuredGrid.h>
#include <vtk-9.0/vtkCellData.h>
#include <vtk-9.0/vtkGenericDataObjectReader.h>
#include <vtk-9.0/vtkGenericDataObjectWriter.h>
template<typename vtk_grid>
int ReadFileVtk(const std::string& name_file_vtk, vtkSmartPointer<vtk_grid>& unstructured_grid) {
	
	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk.c_str());
	reader_vtk->Update();

	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructured_grid = reader_vtk->GetUnstructuredGridOutput();
		unstructured_grid->Modified();
	}
	else
	{
		RETURN_ERRS("Error read file_vtk\n file: %s is not UnstructuredGrid\n", name_file_vtk.c_str());
	}
	
	reader_vtk->GetOutput()->GlobalReleaseDataFlagOn();

	std::cout << "Grid has Cell: " << unstructured_grid->GetNumberOfCells() << '\n';
	return 0;
}

template<typename vtk_grid>
size_t ReadFileVtk(const size_t class_file_vtk, const std::string& name_file_vtk, vtkSmartPointer<vtk_grid>& unstructuredgrid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print/*=false*/) 
{
	ReadFileVtk(name_file_vtk, unstructuredgrid);	

	switch (class_file_vtk) 
	{
	case 0:
		density = NULL;
		absorp_coef = NULL;
		rad_en_loose_rate = NULL;
	case 1:
		density = unstructuredgrid->GetCellData()->GetScalars("alpha");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("alpha");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("Q");
		break;
	case 2:
		density = unstructuredgrid->GetCellData()->GetScalars("density");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("absorp_coef");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("radEnLooseRate");
		break;
	case 5:
		density = unstructuredgrid->GetCellData()->GetScalars("density");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("pressure");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("velocity");  // GetScalars("radEnLooseRate");
		break;
	default:
		RETURN_ERR("Bad type vtk\n");
	}

	if (is_print) 
	{
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfPoints() << " points.\n";
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfCells() << " cells.\n";
		if (class_file_vtk) {
			std::cout << "density_Size: " << density->GetSize() << '\n';
			std::cout << "absorp_coef_Size: " << absorp_coef->GetSize() << '\n';
			std::cout << "Q_Size: " << rad_en_loose_rate->GetSize() << '\n';
		}
	}

	
	return 0;
}

#endif //USE_VTK && !READER_VTK
