#[=======================================================================[.rst:
FindHYPRE
-------

Finds the HYPRE library.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``HYPRE_FOUND``
  True if the system has the HYPRE library.
``HYPRE_INCLUDE_DIRS``
  Include directories needed to use HYPRE.
``HYPRE_LIBRARIES``
  Libraries needed to link to HYPRE.
#]=======================================================================]

find_package(HYPRE QUIET NO_MODULE)
if(HYPRE_FOUND)
    message("HYPRE ${HYPRE_FIND_VERSION} found.")
    if(HYPRE_FIND_COMPONENTS)
        message("HYPRE components found:")
        message("${HYPRE_FIND_COMPONENTS}")
    endif()
    return()
endif()

set(HYPRE_FOUND TRUE)

## Find headers and libraries
find_library(HYPRE_LIBRARIES HYPRE)
find_path(HYPRE_INCLUDE_DIRS HYPRE.h)
if(NOT HYPRE_LIBRARIES OR NOT HYPRE_INCLUDE_DIRS)
    set(HYPRE_FOUND FALSE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE DEFAULT_MSG
                                  HYPRE_FOUND
                                  HYPRE_LIBRARIES
                                  HYPRE_INCLUDE_DIRS)

if (HYPRE_FOUND)
    add_library(HYPRE::HYPRE UNKNOWN IMPORTED)
    set_target_properties(HYPRE::HYPRE PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${HYPRE_LIBRARIES})
endif()

mark_as_advanced(HYPRE_LIBRARIES HYPRE_INCLUDE_DIRS HYPRE_FOUND)
