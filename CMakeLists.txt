#--------------------------------------------------------------------
# Preliminary Settings
#--------------------------------------------------------------------

cmake_minimum_required (VERSION 3.16)

project (OpenSMOKEppInterface)

# --- Set module path in order to use custom CMake modules
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{EIGEN_LIBRARY_PATH}/cmake)
set (CMAKE_COLOR_MAKEFILE ON)
set (VERSION "0.2.0")
set (INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation directory for executables")
set (CMAKE_OSX_DEPLOYMENT_TARGET 12.0)
set (CMAKE_CXX_STANDARD 14)

if( NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()

set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set (CMAKE_CXX_FLAGS_RELEASE "-O2")

#--------------------------------------------------------------------
# Compulsory libraries
#--------------------------------------------------------------------

# --- Boost (Compulsory)
find_package (Boost REQUIRED COMPONENTS date_time filesystem program_options system regex timer chrono)
message (STATUS "Boost_INCLUDE_DIRS    = ${Boost_INCLUDE_DIRS}")
message (STATUS "Boost_LIBRARIES       = ${Boost_LIBRARIES}")

# --- Eigen (Compulsory)
find_package (Eigen3 REQUIRED)
message (STATUS "Eigen_LIBRARIES       = ${EIGEN3_INCLUDE_DIR}")

# --- OpenSmoke++ (Compulsory)
set (OPENSMOKE_INCLUDE_DIR $ENV{OPENSMOKE_LATEST})
message (STATUS "OpenSmoke++           = ${OPENSMOKE_INCLUDE_DIR}")

# --- libconfig++
find_package (CONFIG REQUIRED COMPONENTS libconfig++)
message (STATUS "CONFIG++_INCLUDE_DIR  = ${CONFIG++_INCLUDE_DIR}")
message (STATUS "CONFIG++_LIBRARY      = ${CONFIG++_LIBRARY}")

set(COMPULSORY_LIBRARIES ${Boost_LIBRARIES}
                         ${CONFIG++_LIBRARY})

set(COMPULSORY_INCLUDE ${Boost_INCLUDE_DIRS}
                       ${EIGEN3_INCLUDE_DIR}
                       ${CONFIG++_INCLUDE_DIR}
                       ${OPENSMOKE_INCLUDE_DIR})

#--------------------------------------------------------------------
# Add non-compulsory libraries
#--------------------------------------------------------------------

option (OPENSMOKE_USE_OPENBLAS "Activate OpenBLAS support for OpenSMOKEpp" FALSE)
if (OPENSMOKE_USE_OPENBLAS)
  find_package (OpenBLAS COMPONENTS openblas lapack)
  message (STATUS "OpenBLAS_INCLUDE_DIRS = ${OpenBLAS_INCLUDE_DIR}")
  message (STATUS "BLAS_LIBRARIES        = ${OpenBLAS_LIB}")
  list (APPEND COMPULSORY_INCLUDE ${OpenBLAS_INCLUDE_DIR})
  list (APPEND COMPULSORY_LIBRARIES ${OpenBLAS_LIB})
  add_compile_definitions (OPENSMOKE_USE_OPENBLAS=1)
endif()

option (OPENSMOKE_USE_MKL "Activate MKL/BLAS support for OpenSMOKEpp" FALSE)
if (OPENSMOKE_USE_MKL)
  set (BLA_VENDOR Intel10_64lp)
  find_package (BLAS)
  message (STATUS "MKL_DIR       = ${BLAS_INCLUDE_DIRS}")
  message (STATUS "MKL_LIBRARIES = ${BLAS_LIBRARIES}")
  list (APPEND COMPULSORY_INCLUDE ${BLAS_INCLUDE_DIRS})
  list (APPEND COMPULSORY_LIBRARIES ${BLAS_LIBRARIES})
  add_compile_definitions (OPENSMOKE_USE_MKL=1)
endif()

#--------------------------------------------------------------------
# Include flags
#--------------------------------------------------------------------

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/src
                     ${COMPULSORY_INCLUDE})

#--------------------------------------------------------------------
# Library flags
#--------------------------------------------------------------------

LINK_LIBRARIES (${COMPULSORY_LIBRARIES})

#--------------------------------------------------------------------
# Targets to compile
#--------------------------------------------------------------------

add_library (opensmoke SHARED
    src/OpenSMOKE_Interface.cpp
    src/OpenSMOKE_ODESystem_Interface.h)


#--------------------------------------------------------------------
# Testing
#--------------------------------------------------------------------

enable_testing()
include (tests/CMakeTests.cmake)

