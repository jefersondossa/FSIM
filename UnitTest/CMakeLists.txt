enable_testing()

add_library(test_library STATIC test_main.cpp)
# set_property(TARGET test_library)

target_link_libraries(test_library PUBLIC fsarl_library Catch2::Catch2)

add_unit_test(TestShapeFunction TestShapeFunction.cpp)
add_unit_test(TestGeom TestGeom.cpp)

