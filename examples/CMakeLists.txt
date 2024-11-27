add_executable(cantilever_3d_beam cantilever_3d_beam.cpp)
set_property(TARGET cantilever_3d_beam PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(cantilever_3d_beam PRIVATE fsarl_library Eigen3::Eigen triangle ${FOR_FILES} ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})
