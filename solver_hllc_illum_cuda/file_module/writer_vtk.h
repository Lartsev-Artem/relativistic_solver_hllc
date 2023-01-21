#ifndef WRITER_VTK
#define WRITER_VTK

#include "../prj_config.h"

#ifdef USE_VTK
#include "../global_headers.h"
template<typename str_t>
int WriteVtkGrid(str_t name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid> ugrid, const bool format_bin = true)
{
	vtkSmartPointer<vtkGenericDataObjectWriter> writer = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	if(format_bin) writer->SetFileTypeToBinary();
	writer->SetFileName(name_file_vtk.c_str());
	writer->SetInputData(ugrid);
	writer->Write();
	return 0;
}

#endif // USE_VTK
#endif // !WRITER_VTK