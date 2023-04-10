cmake_minimum_required(VERSION 3.8)

macro(package_add_test TESTNAME)
    find_package(GTest REQUIRED)
    add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME} MFEM::MFEM GTest::gtest_main)
    gtest_discover_tests(${TESTNAME}
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

#for tests using mfem
function(add_mfem_test test_name source_file)
add_mfem_executable(${test_name} ${source_file})
gtest_discover_tests(${test_name}
    WORKING_DIRECTORY ${PROJECT_DIR}
    PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
)
set_target_properties(${test_name} PROPERTIES FOLDER tests)
endfunction()

if(NOT DEFINED CMAKE_CUDA_COMPILER OR ${CMAKE_CUDA_COMPILER} MATCHES "NOTFOUND")
package_add_test(AStyleTests AStyleTests.cpp)
endif()

add_mfem_test(Test Test.cpp)