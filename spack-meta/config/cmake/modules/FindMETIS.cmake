#[=======================================================================[.rst:
FindMETIS
-------

Finds the METIS library.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``METIS_FOUND``
  True if the system has the METIS library.
``METIS_INCLUDE_DIRS``
  Include directories needed to use METIS.
``METIS_LIBRARIES``
  Libraries needed to link to METIS.
#]=======================================================================]

find_package(METIS QUIET NO_MODULE)
if(METIS_FOUND)
    message("METIS ${METIS_FIND_VERSION} found.")
    if(METIS_FIND_COMPONENTS)
        message("METIS components found:")
        message("${METIS_FIND_COMPONENTS}")
    endif()
    return()
endif()

set(METIS_FOUND TRUE)

## Find headers and libraries
find_library(METIS_LIBRARIES metis)
find_path(METIS_INCLUDE_DIRS metis.h)
if(NOT METIS_LIBRARIES OR NOT METIS_INCLUDE_DIRS)
    set(METIS_FOUND FALSE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS DEFAULT_MSG
                                  METIS_FOUND
                                  METIS_LIBRARIES
                                  METIS_INCLUDE_DIRS)

if (METIS_FOUND)
    add_library(METIS::METIS UNKNOWN IMPORTED)
    set_target_properties(METIS::METIS PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${METIS_LIBRARIES})
endif()

mark_as_advanced(METIS_LIBRARIES METIS_INCLUDE_DIRS METIS_FOUND)
