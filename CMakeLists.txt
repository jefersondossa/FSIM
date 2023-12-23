cmake_minimum_required (VERSION 3.11.0)
project (FSArl C CXX Fortran)
enable_language( Fortran CXX )

list (APPEND CMAKE_MODULE_PATH "${FSArl_SOURCE_DIR}/cmake-modules")

find_package (MPI REQUIRED)
find_package (PETSc REQUIRED)   
find_package (METIS REQUIRED)
find_package (HDF5 REQUIRED)
# find_package (Gmsh REQUIRED)

# Essential: include our directories first otherwise we can get internal headers from some installed path
include_directories(include ${MPI_INCLUDE_PATH} ${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})

include_directories ("${FSArl_SOURCE_DIR}/include" "${FSArl_SOURCE_DIR}/include" ${PETSC_INCLUDES} ${PETSC_KSP_LIB})
add_definitions (${PETSC_DEFINITIONS})
set (FSArl_DEPENDENT_LIBRARIES "${PETSC_LIBRARIES}")
set (FSArl_LIBRARIES FSArl ${FSArl_DEPENDENT_LIBRARIES})
include_directories (${HDF5_INCLUDE_DIRS} ${HDF5_DEFINITIONS} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
include_directories(${MPI_INCLUDE_PATH})

include(FortranCInterface)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#[[
If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their definition
and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
to Debug prior to calling PROJECT(). This assumes you're only build one type of solution
a time. If you're working with multiple solutions like Visual Studio Professional can do,
you need to comment this out.
]]
if(DEFINED CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: Debug, Release.")
else() # CMAKE_BUILD_TYPE
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: Debug, Release.")
endif() # CMAKE_BUILD_TYPE
 
set(CMAKE_CXX_STANDARD 17)

# Eigen
include(cmake-modules/EnableEigen.cmake)
enable_eigen()

add_subdirectory(triangle)
add_subdirectory(Util)
add_subdirectory(Analysis)
add_subdirectory(Integral)
add_subdirectory(Mesh)
add_subdirectory(Plasticity)
add_subdirectory(PostProcess)
add_subdirectory(PreProcess)
add_subdirectory(Shape)
add_subdirectory(WeakForm)
# add_subdirectory(examples)
# add_subdirectory(main)

#Specify where the headers files must be found
include_directories(
	"${PROJECT_SOURCE_DIR}/Analysis"
	"${PROJECT_SOURCE_DIR}/Integral"
	"${PROJECT_SOURCE_DIR}/Mesh"
	"${PROJECT_SOURCE_DIR}/Plasticity"
	"${PROJECT_SOURCE_DIR}/PostProcess"
	"${PROJECT_SOURCE_DIR}/PreProcess"
	"${PROJECT_SOURCE_DIR}/Shape"
	"${PROJECT_SOURCE_DIR}/WeakForm"
	"${PROJECT_SOURCE_DIR}/Util"
	"${PROJECT_SOURCE_DIR}/triangle"
	"${PROJECT_SOURCE_DIR}/WeakForm/Arlequin"
	"${PROJECT_SOURCE_DIR}/WeakForm/Elasticity"
	"${PROJECT_SOURCE_DIR}/WeakForm/Poisson"
	"${PROJECT_SOURCE_DIR}/WeakForm/NavierStokes"
	"${PROJECT_SOURCE_DIR}/WeakForm/Stokes"
)

file(GLOB MY_MESH_FILES
  "*.msh" "*.geo" "*.txt"
)
file(COPY ${MY_MESH_FILES} DESTINATION .)

file(GLOB for_files *.for)
set(FOR_FILES ${for_files})

add_library(fsarl_library ${UTIL} ${ANALYSIS} ${INTEGRAL} ${MESH} ${PLASTICITY} ${POSTPROCESS} ${PREPROCESS} ${SHAPE} ${WEAKFORM})
set_property(TARGET fsarl_library PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(fsarl_library PRIVATE Eigen3::Eigen ${FOR_FILES})

add_executable(f main.cpp )
set_property(TARGET f PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries (f fsarl_library Eigen3::Eigen triangle ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})


