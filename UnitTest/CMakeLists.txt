enable_testing()

add_library(test_library STATIC test_main.cpp)
# set_property(TARGET test_library)

target_link_libraries(test_library PUBLIC fsarl_library Catch2::Catch2 Eigen3::Eigen)

add_unit_test(TestShapeFunction TestShapeFunction.cpp)
add_unit_test(TestGeom TestGeom.cpp)
add_unit_test(TestTensor TestTensor.cpp)


add_executable(testTruss testTruss.cpp)
target_link_libraries(testTruss PRIVATE fsarl_library ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})