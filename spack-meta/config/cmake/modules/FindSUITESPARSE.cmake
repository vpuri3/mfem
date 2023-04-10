#[=======================================================================[.rst:
FindSUITESPARSE
-------

Finds the SUITESPARSE library.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``SUITESPARSE_FOUND``
  True if the system has the SUITESPARSE library.
``SUITESPARSE_INCLUDE_DIRS``
  Include directories needed to use SUITESPARSE.
``SUITESPARSE_LIBRARIES``
  Libraries needed to link to SUITESPARSE.
#]=======================================================================]

find_package(SUITESPARSE QUIET NO_MODULE)
if(SUITESPARSE_FOUND)
    message("SUITESPARSE ${SUITESPARSE_FIND_VERSION} found.")
    if(SUITESPARSE_FIND_COMPONENTS)
        message("SUITESPARSE components found:")
        message("${SUITESPARSE_FIND_COMPONENTS}")
    endif()
    return()
endif()

set(SUITESPARSE_FOUND TRUE)

## Find headers and libraries
find_library(SUITESPARSE_CONFIG_LIBRARIES suitesparseconfig)
find_path(SUITESPARSE_CONFIG_INCLUDE_DIR SuiteSparse_config.h)
if(NOT SUITESPARSE_CONFIG_LIBRARIES OR NOT SUITESPARSE_CONFIG_INCLUDE_DIR)
    set(SUITESPARSE_FOUND FALSE)
endif()

find_library(UMFPACK_LIBRARIES umfpack)
find_path(UMFPACK_INCLUDE_DIR umfpack.h)
if(NOT UMFPACK_LIBRARIES OR NOT UMFPACK_INCLUDE_DIR)
    set(SUITESPARSE_FOUND FALSE)
endif()

find_library(KLU_LIBRARIES klu)
find_path(KLU_INCLUDE_DIR klu.h)
if(NOT KLU_LIBRARIES OR NOT KLU_INCLUDE_DIR)
    set(SUITESPARSE_FOUND FALSE)
endif()

list(APPEND SUITESPARSE_LIBRARIES
    ${SUITESPARSE_CONFIG_LIBRARIES}
    ${UMFPACK_LIBRARIES}
    ${KLU_LIBRARIES}
    )

list(APPEND SUITESPARSE_INCLUDE_DIRS
    ${SUITESPARSE_CONFIG_INCLUDE_DIR}
    ${UMFPACK_INCLUDE_DIR}
    ${KLU_INCLUDE_DIR}
    )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUITESPARSE DEFAULT_MSG
                                  SUITESPARSE_FOUND
                                  SUITESPARSE_LIBRARIES
                                  SUITESPARSE_INCLUDE_DIRS)

if (SUITESPARSE_FOUND)
    add_library(SUITESPARSE::SUITESPARSE UNKNOWN IMPORTED)
    set_target_properties(SUITESPARSE::SUITESPARSE PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${SUITESPARSE_CONFIG_INCLUDE_DIR}"
        IMPORTED_LOCATION ${SUITESPARSE_CONFIG_LIBRARIES})
    add_library(SUITESPARSE::UMFPACK UNKNOWN IMPORTED)
    set_target_properties(SUITESPARSE::UMFPACK PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIR}"
        IMPORTED_LOCATION ${UMFPACK_LIBRARIES})
    add_library(SUITESPARSE::KLU UNKNOWN IMPORTED)
    set_target_properties(SUITESPARSE::KLU PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${KLU_INCLUDE_DIR}"
        IMPORTED_LOCATION ${KLU_LIBRARIES})
endif()

mark_as_advanced(SUITESPARSE_LIBRARIES SUITESPARSE_INCLUDE_DIRS SUITESPARSE_FOUND)
