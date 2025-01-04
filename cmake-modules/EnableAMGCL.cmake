# Here we will search for amgcl. If amgcl is not found, we will download it.
function(enable_amgcl)
    #perhaps amgcl was already downloaded when running cmake
    if(NOT amgcl_POPULATED)
        find_package(amgcl QUIET CONFIG)

        if(NOT AMGCL_FOUND)
            # Couldn't load via target, so fall back to allowing module mode finding, which will pick up
            # tools/Findamgcl.cmake
            find_package(amgcl QUIET)
        endif()
        if(NOT AMGCL_FOUND)
            set(AMGCL_VERSION_STRING "1.4.4")
            include(FetchContent)
            FetchContent_Declare(
                amgcl
                GIT_REPOSITORY https://github.com/ddemidov/amgcl.git
                GIT_TAG ${AMGCL_VERSION_STRING})

            FetchContent_GetProperties(amgcl)
            if(NOT amgcl_POPULATED)
                message(STATUS "Downloading amgcl")
                FetchContent_Populate(amgcl)
                set(DOWNLOADED_EIGEN TRUE)
            endif()

            set(AMGCL_INCLUDE_DIR ${amgcl_SOURCE_DIR})
            set(AMGCL_FOUND TRUE)
        endif()
    endif()
    if(AMGCL_FOUND)
        # if amgcl was downloaded, the target was not created
        if(NOT TARGET amgcl::amgcl)
            add_library(amgcl::amgcl IMPORTED INTERFACE)
            set_property(TARGET amgcl::amgcl PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                "${AMGCL_INCLUDE_DIR}")
        endif()

        if(NOT AMGCL_VERSION AND AMGCL_VERSION_STRING)
            set(AMGCL_VERSION ${AMGCL_VERSION_STRING})
        endif()
        message(STATUS "amgcl found at ${AMGCL_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR "Could not satisfy dependency: amgcl")
    endif()
endfunction()