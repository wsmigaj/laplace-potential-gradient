cmake_minimum_required(VERSION 2.8)

# Find the BEM++ library installed on your disk
find_package(Bempp REQUIRED CONFIG)
find_library(BEMPP_LIBRARY bempp PATHS ${BEMPP_LIBRARY_DIR})
# Add the BEM++ installation directory to CMake search path to ensure that
# e.g. the version of SWIG coming with BEM++ is used
set(CMAKE_PREFIX_PATH "${BEMPP_PREFIX}/bempp;${CMAKE_PREFIX_PATH}"
    CACHE STRING "")
# Make the new module be installed in the main BEM++ directory, by default
set(CMAKE_INSTALL_PREFIX "${BEMPP_PREFIX}/bempp" CACHE PATH "" FORCE)
# By default, compile the new module with optimizations turned on
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
       FORCE)
endif ()

# The laplace_gradient library
add_library(laplace_gradient SHARED 
    laplace_3d_single_layer_potential_operator_gradient.cpp
)
target_link_libraries(laplace_gradient ${BEMPP_LIBRARY} ${BEMPP_TEUCHOS_LIBRARY})
set_target_properties(laplace_gradient PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/bempp/lib")
install(TARGETS laplace_gradient LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/bempp/lib")


# Find SWIG
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS "-modern")

# Add include directories necessary for compilation of the new module
include_directories(${BEMPP_PYTHON_INCLUDE_DIRS}) # Python and NumPy headers
include_directories(${BEMPP_INCLUDE_DIR})
include_directories(${BEMPP_INCLUDE_DIR}/bempp)
include_directories(${BEMPP_INCLUDE_DIR}/bempp/swig)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Specify the source files and name of the new module
set_source_files_properties(laplace_gradient.i PROPERTIES CPLUSPLUS ON)
swig_add_module(laplace_gradient python laplace_gradient.i)
swig_link_libraries(laplace_gradient
    laplace_gradient
    ${BEMPP_PYTHON_LIBRARY}
    ${BEMPP_LIBRARY}
    ${BEMPP_TEUCHOS_LIBRARY})

# Specify files to be installed
install(FILES 
laplace_3d_single_layer_potential_gradient_kernel_functor.hpp
laplace_3d_single_layer_potential_operator_gradient.hpp
    DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bempp/assembly")
set(FILES_TO_INSTALL
    "${CMAKE_BINARY_DIR}/laplace_gradient.py"
    "${CMAKE_BINARY_DIR}/_laplace_gradient.so")
install(FILES ${FILES_TO_INSTALL}
    DESTINATION "${CMAKE_INSTALL_PREFIX}/python/bempp")
