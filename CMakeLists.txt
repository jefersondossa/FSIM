cmake_minimum_required (VERSION 3.20.0)
project (FSArl C CXX Fortran)
enable_language( Fortran CXX )

list (APPEND CMAKE_MODULE_PATH "${FSArl_SOURCE_DIR}/cmake-modules")

option(ENABLE_MPI_PETSC "Compile MPI and PETSc packages " OFF)
if(ENABLE_MPI_PETSC)
   # find_package(MPI REQUIRED)
    find_package(PETSc REQUIRED)
endif()

option(ENABLE_GRAPHINTERFACE "Compile Graphic Interface " OFF)
if(ENABLE_GRAPHINTERFACE) 
	find_package(FLTK REQUIRED)
endif()

option(ENABLE_MATPLOTLIBCPP "Compile MatPlotLib for CPP " OFF)
if(ENABLE_MATPLOTLIBCPP) 
	find_package(Python2 COMPONENTS Development)
endif()

#hide warnings
add_definitions(-w) 

option(ENABLE_METIS "Compile METIS package " OFF)
if (ENABLE_METIS)
	find_package (METIS REQUIRED)
endif()

option(ENABLE_AMGCL "Compile AMGCL package " OFF)
if(ENABLE_AMGCL)
	find_package(OpenMP REQUIRED)
	find_package(Boost REQUIRED)

	include(cmake-modules/EnableAMGCL.cmake)
	enable_amgcl()
endif()

option(ENABLE_HDF5 "Compile HDF5 package" OFF)
if (ENABLE_HDF5)
	find_package (HDF5 REQUIRED)
endif()
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

# Catch2
include(cmake-modules/EnableCatch2.cmake)
enable_catch2()

# Eigen
include(cmake-modules/EnableEigen.cmake)
enable_eigen()

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

add_subdirectory(dependencies)

add_subdirectory(triangle)
add_subdirectory(Util)
add_subdirectory(Analysis)
add_subdirectory(Integral)
add_subdirectory(Matrix)
add_subdirectory(Mesh)
add_subdirectory(Plasticity)
add_subdirectory(Plot)
add_subdirectory(PostProcess)
add_subdirectory(PreProcess)
add_subdirectory(Shape)
add_subdirectory(Solver)
add_subdirectory(WeakForm)
add_subdirectory(topology_optimization)
if(ENABLE_GRAPHINTERFACE)
	add_subdirectory(GraphicInterface)
endif()
# add_subdirectory(examples)
# add_subdirectory(main)

#Specify where the headers files must be found
set(INCLUDE_DIRS_LIST 
	"${PROJECT_SOURCE_DIR}/Analysis"
	"${PROJECT_SOURCE_DIR}/GraphicInterface"
	"${PROJECT_SOURCE_DIR}/Integral"
	"${PROJECT_SOURCE_DIR}/Matrix"
	"${PROJECT_SOURCE_DIR}/Mesh"
	"${PROJECT_SOURCE_DIR}/Plasticity"
	"${PROJECT_SOURCE_DIR}/Plot"
	"${PROJECT_SOURCE_DIR}/PostProcess"
	"${PROJECT_SOURCE_DIR}/PreProcess"
	"${PROJECT_SOURCE_DIR}/Shape"
	"${PROJECT_SOURCE_DIR}/Solver"
	"${PROJECT_SOURCE_DIR}/WeakForm"
	"${PROJECT_SOURCE_DIR}/Util"
	"${PROJECT_SOURCE_DIR}/triangle"
	"${PROJECT_SOURCE_DIR}/WeakForm/Arlequin"
	"${PROJECT_SOURCE_DIR}/WeakForm/Elasticity"
	"${PROJECT_SOURCE_DIR}/WeakForm/Poisson"
	"${PROJECT_SOURCE_DIR}/WeakForm/PhaseField"
	"${PROJECT_SOURCE_DIR}/WeakForm/NavierStokes"
	"${PROJECT_SOURCE_DIR}/WeakForm/Stokes"
	"${PROJECT_SOURCE_DIR}/WeakForm/UnitTest")

include_directories(${INCLUDE_DIRS_LIST})

file(GLOB MY_MESH_FILES
  "*.msh" "*.geo" "*.txt"
)
file(COPY ${MY_MESH_FILES} DESTINATION .)

file(GLOB for_files *.for)
set(FOR_FILES ${for_files})

add_library(fsarl_library ${UTIL} ${ANALYSIS} ${GRAPHICINTERFACE} ${INTEGRAL} ${MATRIX} ${MESH} ${PLASTICITY} ${POSTPROCESS} ${PREPROCESS} ${SHAPE} ${SOLVER} ${WEAKFORM})
set_property(TARGET fsarl_library PROPERTY LINKER_LANGUAGE Fortran)

target_link_libraries(fsarl_library PUBLIC Eigen3::Eigen ${FLTK_LIBRARIES})
target_include_directories(fsarl_library PUBLIC ${INCLUDE_DIRS_LIST})

target_compile_definitions(fsarl_library
                           PUBLIC
                           "$<$<CONFIG:RELEASE>:RELEASE_BUILD>"
                           "$<$<CONFIG:DEBUG>:DEBUG_BUILD>")
if(ENABLE_MPI_PETSC)
	target_compile_definitions(fsarl_library
		PUBLIC
		"$<$<CONFIG:RELEASE>:HAS_PETSC>"
		"$<$<CONFIG:DEBUG>:HAS_PETSC>")
endif()

if(ENABLE_HDF5)
	target_compile_definitions(fsarl_library
		PUBLIC
		"$<$<CONFIG:RELEASE>:HAS_HDF5>"
		"$<$<CONFIG:DEBUG>:HAS_HDF5>")
endif()

if(ENABLE_AMGCL)
	target_compile_definitions(fsarl_library
		PUBLIC
		"$<$<CONFIG:RELEASE>:HAS_AMGCL>"
		"$<$<CONFIG:DEBUG>:HAS_AMGCL>")

	target_link_libraries(fsarl_library PRIVATE amgcl::amgcl OpenMP::OpenMP_CXX)
endif()

if(ENABLE_GRAPHINTERFACE)
	target_compile_definitions(fsarl_library
		PUBLIC
		"$<$<CONFIG:RELEASE>:HAS_GRAPHINTERFACE>"
		"$<$<CONFIG:DEBUG>:HAS_GRAPHINTERFACE>")
endif()

if(ENABLE_METIS)
	target_compile_definitions(fsarl_library
		PUBLIC
		"$<$<CONFIG:RELEASE>:HAS_METIS>"
		"$<$<CONFIG:DEBUG>:HAS_METIS>")
	target_link_libraries(fsarl_library PRIVATE metis)
endif()

if (ENABLE_MATPLOTLIBCPP)
	message("\n\n
			**************************************************************\n
			To run MATPLOTLIBCPP please provide CMAKE version minimal 3.29\n
			**************************************************************\n\n\n")
	add_library(plot_library ${PLOT})
	target_include_directories(plot_library PRIVATE ${Python2_INCLUDE_DIRS})
	target_link_libraries(plot_library Python2::Python)

	# add_executable(plotexample quiver.cpp)
	# # set_property(TARGET plotexample PROPERTY )
	# target_link_libraries(plotexample plot_library Python2::Python Python2::NumPy)
endif()

add_executable(f main.cpp )
set_property(TARGET f PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries (f fsarl_library Eigen3::Eigen triangle ${FOR_FILES} ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES})

if(ENABLE_GRAPHINTERFACE)
	add_executable(window Window.cpp )
	set_property(TARGET window PROPERTY LINKER_LANGUAGE Fortran)
	target_link_libraries (window fsarl_library Eigen3::Eigen ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} ${HDF5_LIBRARIES} ${FLTK_LIBRARIES} ${OPENGL_LIBRARIES})
endif()

option(COMPILE_EXAMPLES "Compile examples" OFF)
if (COMPILE_EXAMPLES)
	add_subdirectory(examples)
endif()

#To enable unit test using Catch2 UnitTest framework
option(BUILD_UNITTEST "Whether the Catch2 UnitTest will be linked in" OFF)
if(BUILD_UNITTEST)
  include(cmake-modules/add_unit_test.cmake)
  include(cmake-modules/EnableCatch2.cmake)
  enable_catch2()
endif()


if(BUILD_UNITTEST)
    add_subdirectory(UnitTest)
    #Enabling tests
    ENABLE_TESTING()
    include(CTest)
endif()
