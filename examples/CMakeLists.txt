add_executable(cantilever_3d_beam cantilever_3d_beam.cpp)
target_link_libraries(cantilever_3d_beam PRIVATE fsarl_library triangle ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(windAroundBuilding windAroundBuilding.cpp)
# set_property(TARGET windAroundBuilding PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(windAroundBuilding PRIVATE fsarl_library Eigen3::Eigen triangle ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(cantilever_2d_beam cantilever_2d_beam.cpp)
target_link_libraries(cantilever_2d_beam PRIVATE fsarl_library Eigen3::Eigen triangle ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(eso_cantilever_3d_beam eso_cantilever_3d_beam.cpp)
target_link_libraries(eso_cantilever_3d_beam PRIVATE fsarl_library ESO triangle ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(eso_cantilever_2d_beam eso_cantilever_2d_beam.cpp)
target_link_libraries(eso_cantilever_2d_beam PRIVATE fsarl_library ESO triangle ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(basic_2d_phase_field basic_2d_phase_field.cpp)
target_link_libraries(basic_2d_phase_field PRIVATE fsarl_library ESO triangle ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(basic_2d_phase_field_3D basic_2d_phase_field_3D.cpp)
target_link_libraries(basic_2d_phase_field_3D PRIVATE fsarl_library ESO triangle ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(basic_2d_phase_field_fem basic_2d_phase_field_fem.cpp)
target_link_libraries(basic_2d_phase_field_fem PRIVATE fsarl_library ESO triangle vtu_reader ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(basic_2d_phase_field_fem2 basic_2d_phase_field_fem2.cpp)
target_link_libraries(basic_2d_phase_field_fem2 PRIVATE fsarl_library ESO triangle ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(read_vtu read_vtu.cpp)
target_link_libraries(read_vtu PRIVATE vtu_reader)
