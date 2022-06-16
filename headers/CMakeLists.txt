
add_subdirectory(mesh_interface)

file(GLOB headers *.h)

set(HEADERS ${headers} PARENT_SCOPE)
set(MESH_INT ${mesh_int} PARENT_SCOPE)


