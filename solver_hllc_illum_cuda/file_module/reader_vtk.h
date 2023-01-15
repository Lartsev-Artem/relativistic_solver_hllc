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
	reader_vtk->SetFileName(name_file_vtk);
	reader_vtk->Update();

	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructured_grid = reader_vtk->GetUnstructuredGridOutput();
		unstructured_grid->Modified();
	}
	else
	{
		RETURN_ERR("Error read file_vtk\n file: %s is not UnstructuredGrid\n", name_file_vtk.c_str());
	}

	std::cout << "Grid has Cell: " << unstructured_grid->GetNumberOfCells() << '\n';
	return 0;
}
#endif //USE_VTK && !READER_VTK
