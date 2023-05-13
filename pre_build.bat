echo off
set base_dir=D:\Desktop\FilesCourse\Test\
set utils_dir=D:\MiniProgramm\utils_new\
set np=%1
set name=cone
set build=D:\BMSTU\solver_hllc_illum_cuda\x64\Release\build_cone_jet.exe
set set_build_file=D:\Desktop\FilesCourse\settings_file_build.txt

del %base_dir%%name%_m_epart_%np%.txt

Start /B /wait D:\MiniProgramm\FromNetGenToVTK_volume\Run.exe %base_dir%%name%_netgen.txt %base_dir%%name%.vtk

Start /B /wait %build% %set_build_file%

Start /B /wait %utils_dir%netgen_to_metis.exe %base_dir%%name%_netgen.txt %base_dir%%name%_m.txt
Start /B /wait D:\Programm\Metis\METIS-master\build\programs\Release\mpmetis.exe %base_dir%%name%_m.txt %np%

ren %base_dir%%name%_m.txt.epart.%np% %name%_m_epart_%np%.txt
del %base_dir%%name%_m.txt.npart.%np%

Start /B /wait %utils_dir%make_simpleTxtFromData.exe %base_dir%%name%_m_epart_%np%.txt %base_dir%%name%_m_epart_%np%.txt

Start /B /wait %utils_dir%reorder_metis_nodes.exe %base_dir%%name%_m_epart_%np%.txt %base_dir%centers.bin

Start /B /wait %utils_dir%get_graph_metis_to_renumber.exe %base_dir%%name%_m_epart_%np%.txt %base_dir%%name%_m_graph.bin
del %base_dir%%name%_m_graph.bin.txt

ren %base_dir%%name%.vtk %name%_old.vtk

Start /B /wait %utils_dir%re_numbering_grid.exe %base_dir%%name%_old.vtk %base_dir%%name%.vtk %base_dir%%name%_m_graph.bin
Start /B /wait %build% %set_build_file%

Start /B /wait %utils_dir%second_graph_metis_to_mpi.exe %base_dir%%name%_m_epart_%np%.txt %base_dir%pairs.bin %base_dir%%name%_m_graph.bin

ren %base_dir%%name%.vtk %name%_old1.vtk

Start /B /wait %utils_dir%re_numbering_grid.exe %base_dir%%name%_old1.vtk %base_dir%%name%.vtk %base_dir%%name%_m_graph.bin
Start /B /wait %build% %set_build_file%

Start /B /wait %utils_dir%build_geo_files.exe %base_dir%
Start /B /wait %utils_dir%get_mpi_conf.exe %base_dir%geo_faces.bin %base_dir%%name%_m_epart_%np%.txt %base_dir%mpi_conf.txt

del %base_dir%%name%_old.vtk
del %base_dir%%name%_old1.vtk
