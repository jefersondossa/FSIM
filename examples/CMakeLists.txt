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
