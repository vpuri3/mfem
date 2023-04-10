#[=======================================================================[.rst:
FindCEED
-------

Finds the CEED library.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``CEED_FOUND``
  True if the system has the CEED library.
``CEED_INCLUDE_DIRS``
  Include directories needed to use CEED.
``CEED_LIBRARIES``
  Libraries needed to link to CEED.
#]=======================================================================]

find_package(libCEED QUIET NO_MODULE)
if(CEED_FOUND)
    message("CEED ${CEED_FIND_VERSION} found.")
    if(CEED_FIND_COMPONENTS)
        message("CEED components found:")
        message("${CEED_FIND_COMPONENTS}")
    endif()
    return()
endif()

set(CEED_FOUND TRUE)

## Find headers and libraries
find_library(CEED_LIBRARIES ceed)
find_path(CEED_INCLUDE_DIRS ceed.h)
if(NOT CEED_LIBRARIES OR NOT CEED_INCLUDE_DIRS)
    set(CEED_FOUND FALSE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CEED DEFAULT_MSG
                                  CEED_FOUND
                                  CEED_LIBRARIES
                                  CEED_INCLUDE_DIRS)

if (CEED_FOUND)
    add_library(CEED::CEED UNKNOWN IMPORTED)
    set_target_properties(CEED::CEED PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CEED_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${CEED_LIBRARIES}
        )
endif()

mark_as_advanced(CEED_LIBRARIES CEED_INCLUDE_DIRS CEED_FOUND)
