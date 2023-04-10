#[=======================================================================[.rst:
FindMFEM
-------

Finds the MFEM library.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``MFEM_FOUND``
  True if the system has the MFEM library.
``MFEM_INCLUDE_DIRS``
  Include directories needed to use MFEM.
``MFEM_LIBRARIES``
  Libraries needed to link to MFEM.
#]=======================================================================]

find_package(MFEM QUIET NO_MODULE)
if(MFEM_FOUND)
    message("MFEM ${MFEM_FIND_VERSION} found.")
    if(MFEM_FIND_COMPONENTS)
        message("MFEM components found:")
        message("${MFEM_FIND_COMPONENTS}")
    endif()
    return()
endif()

set(MFEM_FOUND TRUE)

## Find headers and libraries
find_library(MFEM_LIBRARIES mfem)
find_path(MFEM_INCLUDE_DIRS mfem.hpp)
if(NOT MFEM_LIBRARIES OR NOT MFEM_INCLUDE_DIRS)
    set(MFEM_FOUND FALSE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MFEM DEFAULT_MSG
                                  MFEM_FOUND
                                  MFEM_LIBRARIES
                                  MFEM_INCLUDE_DIRS)

if (MFEM_FOUND)
    add_library(MFEM::MFEM UNKNOWN IMPORTED)
    set_target_properties(MFEM::MFEM PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${MFEM_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${MFEM_LIBRARIES}
        )
endif()

mark_as_advanced(MFEM_LIBRARIES MFEM_INCLUDE_DIRS MFEM_FOUND)
