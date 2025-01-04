add_library(ESO eso.cpp)
target_include_directories(ESO SYSTEM PUBLIC .)
target_link_libraries(ESO PRIVATE fsarl_library)
